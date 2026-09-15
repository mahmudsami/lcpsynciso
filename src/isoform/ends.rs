//! Step 3: split a structure group into isoforms by where its reads start and end.
//!
//! Every read gets a start key and an end key. With `end_modes`, the keys are peaks found in
//! the reads' start (or end) coordinates: peaks are found with the narrow `peak_width`,
//! because real boundaries can sit ~20 bp apart, and a read joins the nearest peak within
//! the wider `boundary_tol`. Without `end_modes`, the keys are cells of a fixed
//! `boundary_tol` grid.
//!
//! A (start, end) key pair held by at least `min_iso` reads is an isoform. Any other read,
//! typically a truncated copy, joins the best-supported isoform whose extent contains it,
//! or is dropped.

use std::collections::HashMap;

use super::{stats, Cfg, Member};

/// Member read indices of each isoform found in `members`, in no particular order.
pub(super) fn split_by_ends(members: &[Member], cfg: &Cfg) -> Vec<Vec<usize>> {
    let tol = cfg.boundary_tol.max(1) as i32;

    let (start_key, end_key): (Vec<Option<i32>>, Vec<Option<i32>>) = if cfg.end_modes {
        let width = cfg.peak_width.max(1) as i32;
        let start_peaks = find_peaks(members.iter().map(|m| m.1), width, cfg.min_iso);
        let end_peaks = find_peaks(members.iter().map(|m| m.2), width, cfg.min_iso);
        (
            members.iter().map(|m| nearest_peak(&start_peaks, m.1, tol)).collect(),
            members.iter().map(|m| nearest_peak(&end_peaks, m.2, tol)).collect(),
        )
    } else {
        (
            members.iter().map(|m| Some(m.1.div_euclid(tol))).collect(),
            members.iter().map(|m| Some(m.2.div_euclid(tol))).collect(),
        )
    };
    let key_of = |mi: usize| match (start_key[mi], end_key[mi]) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => None,
    };

    let mut cells: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for mi in 0..members.len() {
        if let Some(key) = key_of(mi) {
            cells.entry(key).or_default().push(mi);
        }
    }

    // Dense cells are isoforms, each described by its extent and read count. The extent is
    // the q-th smallest start and q-th largest end, not the median, so that the containment
    // test below admits every read that belongs. q stays well below the cell size; at
    // q = min_iso a minimal cell would shrink to its members' intersection.
    let mut isoforms: HashMap<(i32, i32), (i32, i32, usize)> = HashMap::new();
    for (key, mis) in &cells {
        if mis.len() >= cfg.min_iso {
            let mut starts: Vec<i32> = mis.iter().map(|&mi| members[mi].1).collect();
            starts.sort_unstable();
            let mut ends: Vec<i32> = mis.iter().map(|&mi| members[mi].2).collect();
            ends.sort_unstable();
            let n = mis.len();
            let q = (cfg.min_iso / 2).max(1).min(n);
            isoforms.insert(*key, (starts[q - 1], ends[n - q], n));
        }
    }
    if isoforms.is_empty() {
        // No dense cell (a smooth smear of ends): the whole group is one isoform.
        return vec![members.iter().map(|&(ri, _, _)| ri).collect()];
    }

    let mut out: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for (mi, &(ri, s, e)) in members.iter().enumerate() {
        if let Some(key) = key_of(mi) {
            if isoforms.contains_key(&key) {
                stats::inc(&stats::N_DENSE);
                out.entry(key).or_default().push(ri);
                continue;
            }
        }
        // Otherwise: the best-supported isoform whose extent contains [s, e], within tol.
        let mut best: Option<(i32, i32)> = None;
        let mut best_sup = 0usize;
        for (key, &(is, ie, sup)) in &isoforms {
            if is <= s + tol && ie >= e - tol && sup > best_sup {
                best_sup = sup;
                best = Some(*key);
            }
        }
        if let Some(key) = best {
            stats::inc(&stats::N_ENCL);
            let (is, ie, _) = isoforms[&key];
            if (e - s) > (ie - is) + tol {
                stats::inc(&stats::N_FOLD_WIDER);
            }
            out.entry(key).or_default().push(ri);
        } else {
            stats::inc(&stats::N_DROP);
        }
    }
    out.into_values().collect()
}

/// Peaks in a set of coordinates, ascending: repeatedly take the `width`-wide window holding
/// the most remaining points, record its median, and remove those points. Stops once the
/// densest window has fewer than `min_support` points, since what remains is a smear.
fn find_peaks(xs: impl Iterator<Item = i32>, width: i32, min_support: usize) -> Vec<i32> {
    let mut pts: Vec<i32> = xs.collect();
    pts.sort_unstable();

    let mut peaks: Vec<i32> = Vec::new();
    while !pts.is_empty() {
        // Densest window pts[lo..=hi], by a two-pointer scan.
        let (mut lo, mut hi, mut best) = (0usize, 0usize, 0usize);
        let mut j = 0usize;
        for i in 0..pts.len() {
            if j < i {
                j = i;
            }
            while j + 1 < pts.len() && pts[j + 1] - pts[i] <= width {
                j += 1;
            }
            if j - i + 1 > best {
                best = j - i + 1;
                lo = i;
                hi = j;
            }
        }
        if best < min_support {
            break;
        }
        peaks.push(pts[(lo + hi) / 2]);
        pts.drain(lo..=hi);
    }
    peaks.sort_unstable();
    peaks
}

/// The peak nearest to `x` within `radius` (the lower one on a tie), if any.
fn nearest_peak(peaks: &[i32], x: i32, radius: i32) -> Option<i32> {
    peaks
        .iter()
        .copied()
        .filter(|&m| (m - x).abs() <= radius)
        .min_by_key(|&m| ((m - x).abs(), m))
}
