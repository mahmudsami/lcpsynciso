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
//!
//! With `split_starts` off, every read has the same start key, so only 3' ends split a group,
//! and containment is tested at the 3' end only: a read that reaches further 5' than an
//! isoform is a more complete copy of it, not a different transcript. In cDNA most reads are
//! 5'-truncated, so the truncated majority forms tight start peaks while the few full-length
//! reads, whose starts smear over ~100 bp, form none; splitting by starts then made the
//! truncated majority an isoform of its own and dropped the full-length reads.

use std::collections::HashMap;

use super::{stats, Cfg, Member};

/// Member read indices of each isoform found in `members`, in no particular order.
pub(super) fn split_by_ends(members: &[Member], cfg: &Cfg) -> Vec<Vec<usize>> {
    let tol = cfg.boundary_tol.max(1) as i32;

    let (mut start_key, end_key): (Vec<Option<i32>>, Vec<Option<i32>>) = if cfg.end_modes {
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
    if !cfg.split_starts {
        start_key.fill(Some(0));
    }
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
        // Otherwise: the best-supported isoform whose extent contains [s, e], within tol —
        // at the 3' end only when starts don't split.
        let mut best: Option<(i32, i32)> = None;
        let mut best_sup = 0usize;
        for (key, &(is, ie, sup)) in &isoforms {
            let starts_inside = !cfg.split_starts || is <= s + tol;
            if starts_inside && ie >= e - tol && sup > best_sup {
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

#[cfg(test)]
mod tests {
    use super::*;

    /// Sorted member counts of the isoforms `split_by_ends` returns.
    fn sizes(members: &[Member], cfg: &Cfg) -> Vec<usize> {
        let mut s: Vec<usize> = split_by_ends(members, cfg).iter().map(|m| m.len()).collect();
        s.sort_unstable();
        s
    }

    /// Three full-length reads whose starts smear over ~100 bp, and eight 5'-truncated reads
    /// that all start near 800; every read ends at the same polyA site near 2000.
    fn truncated_majority() -> Vec<Member> {
        let mut m: Vec<Member> = vec![(0, 0, 2000), (1, 40, 2003), (2, 95, 1998)];
        for (i, s) in [800, 802, 805, 801, 803, 799, 804, 806].into_iter().enumerate() {
            m.push((3 + i, s, 2000 + i as i32 % 3));
        }
        m
    }

    #[test]
    fn start_split_drops_full_length_reads_of_a_truncated_majority() {
        // The truncated reads form a start peak, the full-length reads none; not contained
        // in the truncated isoform's extent, they are dropped.
        assert_eq!(sizes(&truncated_majority(), &Cfg::default()), vec![8]);
    }

    #[test]
    fn without_start_split_full_length_reads_join_their_isoform() {
        let cfg = Cfg { split_starts: false, ..Cfg::default() };
        assert_eq!(sizes(&truncated_majority(), &cfg), vec![11]);
    }

    #[test]
    fn without_start_split_alternative_polya_sites_still_split() {
        let cfg = Cfg { split_starts: false, ..Cfg::default() };
        let m: Vec<Member> = vec![
            (0, 0, 2000), (1, 300, 2002), (2, 700, 1999), (3, 900, 2001),
            (4, 5, 2600), (5, 350, 2604), (6, 650, 2598), (7, 1000, 2601),
        ];
        assert_eq!(sizes(&m, &cfg), vec![4, 4]);
    }
}
