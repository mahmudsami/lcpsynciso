//! Step 2: split a structure group into isoforms by where its reads start and end.
//!
//! Every read gets a start key and an end key: the peaks found in the reads' start (or end)
//! coordinates. Peaks are found with the narrow `peak_width`, because real boundaries can
//! sit ~20 bp apart, and a read joins the nearest peak within the wider `boundary_tol`.
//!
//! A (start, end) key pair held by at least `min_iso` reads is an isoform. Any other read,
//! typically a truncated copy, joins the best-supported isoform whose extent contains it,
//! or is dropped.
//!
//! How starts count depends on `start_split`. In cDNA most reads are 5'-truncated, so the
//! truncated majority forms tight start peaks while the few full-length reads, whose starts
//! smear over ~100 bp, form none; splitting by every start peak (`On`) made the truncated
//! majority an isoform of its own and dropped the full-length reads.
//!
//!   * `Off`: every read has the same start key, so only 3' ends split a group, and
//!     containment is tested at the 3' end only: a read that reaches further 5' than an
//!     isoform is a more complete copy of it, not a different transcript.
//!   * `Auto`: as `Off`, except that a downstream start peak that looks like a real
//!     transcription start ([`find_validated_starts`]) keys its reads into an isoform of
//!     their own. The full-length isoform is kept either way, so a wrong call costs an extra
//!     isoform, never the complete sequence.

use std::collections::HashMap;

use super::{stats, PlacedRead, Read, IsoformOptions, StartSplit};

/// Start key of the reads at no validated start: the full-length isoform and its truncated
/// copies.
const NO_START: i32 = i32::MIN;
/// A read starting within this many bp of a validated start starts there. Real starts spread
/// over ~20-50 bp, and this was the best peak window in calibration.
const START_RADIUS: i32 = 25;
/// Sharpness compares the starts within `START_RADIUS` with those within this many bp.
const SHARP_WINDOW: i32 = 300;
/// Without a cap signal, a peak this sharp is a real start. SIRV: truncation peaks mostly
/// below 0.3-0.5, real starts that --start-split off merged away at 0.95-1.00.
const MIN_SHARPNESS: f64 = 0.45;
/// With a cap signal, a peak whose reads carry an untemplated 5' G at least this often is a
/// real start. UHRR HiFi: annotated starts ~97-100%, truncation peaks ~50%.
const MIN_CAPPED: f64 = 0.7;
/// A group shows the cap signal when at least this share of its informative reads carry an
/// untemplated G, out of at least `MIN_CAP_READS`. Libraries without template switching
/// (direct RNA, uncapped spike-ins) show almost none.
const CAP_LIBRARY: f64 = 0.2;
const MIN_CAP_READS: usize = 10;

/// Read indices of each isoform found in `group`, and whether it starts at a validated
/// downstream start, in no particular order. `group` sits on `reads[backbone_idx]`.
pub(super) fn split_by_ends(
    group: &[PlacedRead],
    reads: &[Read],
    backbone_idx: usize,
    options: &IsoformOptions,
) -> Vec<(Vec<usize>, bool)> {
    let tol = options.boundary_tol.max(1) as i32;

    let width = options.peak_width.max(1) as i32;
    let start_peaks = find_peaks(group.iter().map(|&(_, start, _)| start), width, options.min_iso);
    let end_peaks = find_peaks(group.iter().map(|&(_, _, end)| end), width, options.min_iso);
    let mut start_key: Vec<Option<i32>> =
        group.iter().map(|&(_, start, _)| nearest_peak(&start_peaks, start, tol)).collect();
    let end_key: Vec<Option<i32>> =
        group.iter().map(|&(_, _, end)| nearest_peak(&end_peaks, end, tol)).collect();
    match options.start_split {
        StartSplit::On => {}
        StartSplit::Off => start_key.fill(Some(NO_START)),
        StartSplit::Auto => {
            let alt_starts = find_validated_starts(group, reads, backbone_idx, options);
            for (key, &(_, start, _)) in start_key.iter_mut().zip(group) {
                *key = Some(nearest_peak(&alt_starts, start, START_RADIUS).unwrap_or(NO_START));
            }
        }
    }
    let key_of = |group_pos: usize| match (start_key[group_pos], end_key[group_pos]) {
        (Some(a), Some(b)) => Some((a, b)),
        _ => None,
    };

    let mut cells: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for group_pos in 0..group.len() {
        if let Some(key) = key_of(group_pos) {
            cells.entry(key).or_default().push(group_pos);
        }
    }

    // Dense cells are isoforms, each described by its extent and read count. The extent is
    // the q-th smallest start and q-th largest end, not the median, so that the containment
    // test below admits every read that belongs. q stays well below the cell size; at
    // q = min_iso a minimal cell would shrink to its reads' intersection.
    let mut isoform_extents: HashMap<(i32, i32), (i32, i32, usize)> = HashMap::new();
    for (key, group_positions) in &cells {
        if group_positions.len() >= options.min_iso {
            let mut starts: Vec<i32> =
                group_positions.iter().map(|&group_pos| group[group_pos].1).collect();
            starts.sort_unstable();
            let mut ends: Vec<i32> =
                group_positions.iter().map(|&group_pos| group[group_pos].2).collect();
            ends.sort_unstable();
            let n = group_positions.len();
            let q = (options.min_iso / 2).max(1).min(n);
            isoform_extents.insert(*key, (starts[q - 1], ends[n - q], n));
        }
    }
    if isoform_extents.is_empty() {
        // No dense cell (a smooth smear of ends): the whole group is one isoform.
        return vec![(group.iter().map(|&(read_idx, _, _)| read_idx).collect(), false)];
    }

    let mut reads_by_isoform: HashMap<(i32, i32), Vec<usize>> = HashMap::new();
    for (group_pos, &(read_idx, start, end)) in group.iter().enumerate() {
        if let Some(key) = key_of(group_pos) {
            if isoform_extents.contains_key(&key) {
                stats::inc(&stats::N_DENSE);
                reads_by_isoform.entry(key).or_default().push(read_idx);
                continue;
            }
        }
        // Otherwise: the best-supported isoform whose extent contains [start, end], within
        // tol. Unless every start peak splits, the 5' side only binds an isoform at a
        // validated start, which takes no read starting upstream of it.
        let mut best: Option<(i32, i32)> = None;
        let mut best_support = 0usize;
        for (key, &(isoform_start, isoform_end, support)) in &isoform_extents {
            let starts_inside = match options.start_split {
                StartSplit::On => isoform_start <= start + tol,
                _ => key.0 == NO_START || isoform_start <= start + START_RADIUS,
            };
            if starts_inside && isoform_end >= end - tol && support > best_support {
                best_support = support;
                best = Some(*key);
            }
        }
        if let Some(key) = best {
            stats::inc(&stats::N_ENCLOSED);
            let (isoform_start, isoform_end, _) = isoform_extents[&key];
            if (end - start) > (isoform_end - isoform_start) + tol {
                stats::inc(&stats::N_ENCLOSED_WIDER);
            }
            reads_by_isoform.entry(key).or_default().push(read_idx);
        } else {
            stats::inc(&stats::N_DROP);
        }
    }
    let auto = options.start_split == StartSplit::Auto;
    reads_by_isoform
        .into_iter()
        .map(|(key, isoform)| (isoform, auto && key.0 != NO_START))
        .collect()
}

/// Start peaks inside the backbone that look like real transcription starts rather than 5'
/// truncation, ascending. The backbone reaches further 5', so each is a candidate alternative
/// start of a longer transcript. A candidate needs `min_iso` reads within `START_RADIUS`, and
///   * if the group shows a cap signal (template switching adds an untemplated G at a capped
///     5' end), `MIN_CAPPED` of them carrying an untemplated G;
///   * otherwise (direct RNA, uncapped spike-ins), a sharp peak: `MIN_SHARPNESS` of the starts
///     within `SHARP_WINDOW`.
/// Thresholds are from UHRR HiFi against GENCODE and SIRV spike-ins against their annotation.
fn find_validated_starts(
    group: &[PlacedRead],
    reads: &[Read],
    backbone_idx: usize,
    options: &IsoformOptions,
) -> Vec<i32> {
    let width = options.peak_width.max(1) as i32;
    let candidates: Vec<i32> =
        find_peaks(group.iter().map(|&(_, start, _)| start), width, options.min_iso)
            .into_iter()
            .filter(|&peak| peak > START_RADIUS)
            .collect();
    if candidates.is_empty() {
        return candidates;
    }
    let backbone_seq = &reads[backbone_idx].seq;
    let capped: Vec<Option<bool>> = group
        .iter()
        .map(|&(read_idx, start, _)| {
            if start > START_RADIUS {
                has_untemplated_g(&reads[read_idx].seq, backbone_seq, start as usize)
            } else {
                None
            }
        })
        .collect();
    let tally = |group_positions: &mut dyn Iterator<Item = usize>| {
        group_positions
            .filter_map(|group_pos| capped[group_pos])
            .fold((0usize, 0usize), |(n, g), c| (n + 1, g + c as usize))
    };
    let (known, with_g) = tally(&mut (0..group.len()));
    let cap_signal = known >= MIN_CAP_READS && with_g as f64 >= CAP_LIBRARY * known as f64;

    let mut validated = Vec::new();
    for peak in candidates {
        let peak_reads: Vec<usize> = (0..group.len())
            .filter(|&group_pos| (group[group_pos].1 - peak).abs() <= START_RADIUS)
            .collect();
        if peak_reads.len() < options.min_iso {
            continue;
        }
        let real = if cap_signal {
            let (n, g) = tally(&mut peak_reads.iter().copied());
            n >= options.min_iso && g as f64 >= MIN_CAPPED * n as f64
        } else {
            let near = group
                .iter()
                .filter(|&&(_, start, _)| (start - peak).abs() <= SHARP_WINDOW)
                .count();
            peak_reads.len() as f64 >= MIN_SHARPNESS * near as f64
        };
        if real {
            validated.push(peak);
        }
    }
    validated
}

/// Whether `read_seq`, whose first matched base sits at backbone position `start`, begins
/// with an untemplated G: 1-3 leading bases, the first a G, that the backbone doesn't have
/// there. None if the read's start doesn't match the backbone exactly nearby (sequencing
/// errors), so the read tells nothing. A G that the transcript itself has just before the
/// start is indistinguishable and counts as templated.
fn has_untemplated_g(read_seq: &[u8], backbone_seq: &[u8], start: usize) -> Option<bool> {
    const L: usize = 12;
    for lead in 0..=3usize {
        for b in [start, start.saturating_sub(1), start + 1] {
            if b + L <= backbone_seq.len()
                && lead + L <= read_seq.len()
                && read_seq[lead..lead + L] == backbone_seq[b..b + L]
            {
                return Some(lead > 0 && read_seq[0] == b'G');
            }
        }
    }
    None
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

    /// Sorted read counts of the isoforms `split_by_ends` returns; `On` and `Off` never
    /// look at read sequences.
    fn sizes(group: &[PlacedRead], options: &IsoformOptions) -> Vec<usize> {
        let reads: Vec<Read> = (0..group.len())
            .map(|read_idx| Read { name: read_idx.to_string(), seq: Vec::new() })
            .collect();
        let mut s: Vec<usize> = split_by_ends(group, &reads, 0, options)
            .iter()
            .map(|isoform| isoform.0.len())
            .collect();
        s.sort_unstable();
        s
    }

    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    /// A 2 kb transcript: three full-length reads (the first is the backbone), six reads
    /// starting at 800, and truncated reads starting every 120 bp from 230, all ending at the
    /// 3' end. `g_at_800` / `g_truncated`: give those reads an untemplated 5' G. `smear`:
    /// add reads starting around 800, so that peak is not sharp.
    fn alt_start_group(g_at_800: bool, g_truncated: bool, smear: bool) -> Vec<(Vec<usize>, bool)> {
        let mut truth = random_seq(2000, 17);
        let mut starts: Vec<usize> = vec![0, 0, 0];
        starts.extend([800; 6]);
        starts.extend((230..1400).step_by(120));
        if smear {
            starts.extend((540..1080).step_by(20).filter(|p| (*p as i32 - 800).abs() > 30));
        }
        for &p in &starts {
            if p > 0 && truth[p - 1] == b'G' {
                truth[p - 1] = b'A'; // so a leading G is untemplated, not the transcript's
            }
        }
        let reads: Vec<Read> = starts
            .iter()
            .enumerate()
            .map(|(read_idx, &p)| {
                let g = read_idx >= 3 && if p == 800 { g_at_800 } else { g_truncated };
                let mut seq = if g { vec![b'G'] } else { Vec::new() };
                seq.extend_from_slice(&truth[p..]);
                Read { name: read_idx.to_string(), seq }
            })
            .collect();
        let group: Vec<PlacedRead> =
            starts.iter().enumerate().map(|(read_idx, &p)| (read_idx, p as i32, 2000)).collect();
        let options = IsoformOptions { start_split: StartSplit::Auto, ..IsoformOptions::default() };
        split_by_ends(&group, &reads, 0, &options)
    }

    #[test]
    fn auto_keeps_a_capped_downstream_start_beside_the_full_length_isoform() {
        let isoforms = alt_start_group(true, true, false);
        let alt: Vec<&Vec<usize>> = isoforms.iter().filter(|i| i.1).map(|i| &i.0).collect();
        assert_eq!(alt.len(), 1);
        assert!((3..9).all(|r| alt[0].contains(&r)), "the six reads at 800 start there");
        assert!(isoforms.iter().any(|i| !i.1 && (0..3).all(|r| i.0.contains(&r))));
    }

    #[test]
    fn auto_treats_an_uncapped_peak_in_a_capped_library_as_truncation() {
        let isoforms = alt_start_group(false, true, false);
        assert!(isoforms.iter().all(|i| !i.1));
        assert_eq!(isoforms.len(), 1);
    }

    #[test]
    fn auto_without_a_cap_signal_keeps_a_sharp_peak_and_folds_a_smeared_one() {
        assert!(alt_start_group(false, false, false).iter().any(|i| i.1));
        assert!(alt_start_group(false, false, true).iter().all(|i| !i.1));
    }

    /// Three full-length reads whose starts smear over ~100 bp, and eight 5'-truncated reads
    /// that all start near 800; every read ends at the same polyA site near 2000.
    fn truncated_majority() -> Vec<PlacedRead> {
        let mut group: Vec<PlacedRead> = vec![(0, 0, 2000), (1, 40, 2003), (2, 95, 1998)];
        for (i, start) in [800, 802, 805, 801, 803, 799, 804, 806].into_iter().enumerate() {
            group.push((3 + i, start, 2000 + i as i32 % 3));
        }
        group
    }

    #[test]
    fn start_split_drops_full_length_reads_of_a_truncated_majority() {
        // The truncated reads form a start peak, the full-length reads none; not contained
        // in the truncated isoform's extent, they are dropped.
        assert_eq!(sizes(&truncated_majority(), &IsoformOptions::default()), vec![8]);
    }

    #[test]
    fn without_start_split_full_length_reads_join_their_isoform() {
        let options = IsoformOptions { start_split: StartSplit::Off, ..IsoformOptions::default() };
        assert_eq!(sizes(&truncated_majority(), &options), vec![11]);
    }

    #[test]
    fn without_start_split_alternative_polya_sites_still_split() {
        let options = IsoformOptions { start_split: StartSplit::Off, ..IsoformOptions::default() };
        let group: Vec<PlacedRead> = vec![
            (0, 0, 2000), (1, 300, 2002), (2, 700, 1999), (3, 900, 2001),
            (4, 5, 2600), (5, 350, 2604), (6, 650, 2598), (7, 1000, 2601),
        ];
        assert_eq!(sizes(&group, &options), vec![4, 4]);
    }
}
