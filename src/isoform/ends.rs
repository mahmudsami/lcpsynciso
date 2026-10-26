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

use std::collections::BTreeMap;

use super::{stats, CapSignal, IsoformOptions, PlacedRead, Read, StartLibrary, StartSplit};

/// Start key of the reads at no validated start: the full-length isoform and its truncated
/// copies.
const NO_START: i32 = i32::MIN;
/// A read starting within this many bp of a validated start starts there. Real starts spread
/// over ~20-50 bp, and this was the best peak window in calibration.
const START_RADIUS: i32 = 25;
/// A validated start is the centre of its start site (a pile's median), so `START_RADIUS`
/// around it is the site. A start that isn't validated is its backbone's start, the most
/// upstream read of its site: the site's edge, which reaches up to the full site width
/// downstream (refTSS: 94% of regions are at most 50 bp wide). A start candidate must lie
/// further than this into the backbone, and merging keeps such a start and a validated one
/// apart only if they differ by more than this.
pub(super) const SITE_WIDTH: i32 = 2 * START_RADIUS;
/// Sharpness compares the starts within `START_RADIUS` with those within this many bp.
const SHARP_WINDOW: i32 = 300;
/// Without a cap signal and without a threshold from the library (`find-isoforms`, or too few
/// start candidates), a peak this sharp is a real start. The library's Otsu cut fell at
/// 0.41-0.49 on SIRV (HiFi, ONT cDNA, ONT direct RNA) and UHRR (same three) start peaks.
const MIN_SHARPNESS: f64 = 0.45;
/// With a cap signal and no threshold from the library (`find-isoforms`, or too few start
/// candidates), a peak whose reads carry an untemplated 5' G at least this often is a real start.
const MIN_CAPPED: f64 = 0.7;
/// A library (with `find-isoforms`, a structure group) shows the cap signal when at least
/// this share of its labelled reads carry an untemplated G, out of at least `MIN_CAP_READS`.
/// Capped share of labelled reads measured with --dump-starts: UHRR PacBio HiFi 41% and ONT
/// cDNA 27% (template and strand switching), ONT direct RNA 12% and SIRV spike-ins 4% (none).
const CAP_LIBRARY: f64 = 0.2;
const MIN_CAP_READS: usize = 10;
/// Outer windows (bp) recorded by --dump-starts, for choosing the sharpness window.
const DUMP_WINDOWS: [i32; 9] = [50, 100, 150, 200, 300, 500, 750, 1000, 1500];
/// Fewest start candidates with labelled reads for setting the library's threshold from
/// their capped shares; with fewer, `MIN_CAPPED` is used.
const MIN_LIBRARY_CANDIDATES: usize = 100;

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

    // Ordered maps throughout: which isoform a tied read joins, and the order isoforms reach
    // merging, then do not depend on hashing, so the output is the same on every run.
    let mut cells: BTreeMap<(i32, i32), Vec<usize>> = BTreeMap::new();
    for group_pos in 0..group.len() {
        if let Some(key) = key_of(group_pos) {
            cells.entry(key).or_default().push(group_pos);
        }
    }

    // Dense cells are isoforms, each described by its extent and read count. The extent is
    // the q-th smallest start and q-th largest end, not the median, so that the containment
    // test below admits every read that belongs. q stays well below the cell size; at
    // q = min_iso a minimal cell would shrink to its reads' intersection.
    let mut isoform_extents: BTreeMap<(i32, i32), (i32, i32, usize)> = BTreeMap::new();
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

    let mut reads_by_isoform: BTreeMap<(i32, i32), Vec<usize>> = BTreeMap::new();
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
/// truncation, ascending. The backbone reaches further 5', so each peak more than
/// `SITE_WIDTH` into it is a candidate alternative start of a longer transcript. A candidate
/// needs `min_iso` reads within `START_RADIUS`, and
///   * if the library shows a cap signal (template switching adds an untemplated G at a capped
///     5' end), the library's share of them carrying an untemplated G (`StartLibrary`);
///   * otherwise (direct RNA, uncapped spike-ins), a sharp peak: `MIN_SHARPNESS` of the starts
///     within `SHARP_WINDOW`.
/// Without a library decision (`find-isoforms`), the cap signal is decided per structure group.
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
            .filter(|&peak| peak > SITE_WIDTH)
            .collect();
    // With --dump-starts every group is labelled, so the dump covers groups without a
    // candidate too; the validated starts are the same.
    let dump = stats::start_dump_on();
    if candidates.is_empty() && !dump {
        return candidates;
    }
    let backbone_seq = &reads[backbone_idx].seq;
    let capped: Vec<Option<bool>> = group
        .iter()
        .map(|&(read_idx, start, _)| {
            if start > START_RADIUS {
                let lead = untemplated_lead(&reads[read_idx].seq, backbone_seq, start as usize);
                if let (true, Some((n, first))) = (dump, lead) {
                    stats::record_lead(n, first);
                }
                lead.map(|(n, first)| n > 0 && first == b'G')
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
    let (cap_signal, min_capped, min_sharpness) = match options.start_library {
        Some(library) => (library.cap_signal, library.threshold, library.sharpness),
        None => (
            match options.cap_signal {
                CapSignal::On => true,
                CapSignal::Off => false,
                CapSignal::Auto => {
                    known >= MIN_CAP_READS && with_g as f64 >= CAP_LIBRARY * known as f64
                }
            },
            options.cap_threshold.unwrap_or(MIN_CAPPED),
            options.sharp_threshold.unwrap_or(MIN_SHARPNESS),
        ),
    };
    if dump {
        stats::start_dump(format!("group\t{}\t{known}\t{with_g}\t{}", group.len(), cap_signal as u8));
    }

    let mut validated = Vec::new();
    for peak in candidates {
        let peak_reads: Vec<usize> = (0..group.len())
            .filter(|&group_pos| (group[group_pos].1 - peak).abs() <= START_RADIUS)
            .collect();
        if peak_reads.len() < options.min_iso {
            continue;
        }
        let (n, g) = tally(&mut peak_reads.iter().copied());
        let near = group
            .iter()
            .filter(|&&(_, start, _)| (start - peak).abs() <= SHARP_WINDOW)
            .count();
        let real = if cap_signal {
            n >= options.min_iso && g as f64 >= min_capped * n as f64
        } else {
            peak_reads.len() as f64 >= min_sharpness * near as f64
        };
        if dump {
            // Reads starting within each of several windows, for choosing the sharpness
            // window; and the backbone, for placing the candidate on an annotation.
            let windows: Vec<String> = DUMP_WINDOWS
                .iter()
                .map(|&w| {
                    let c = group.iter().filter(|&&(_, start, _)| (start - peak).abs() <= w).count();
                    format!("{w}:{c}")
                })
                .collect();
            stats::start_dump(format!(
                "cand\t{peak}\t{}\t{n}\t{g}\t{near}\t{}\t{}\t{}\t{}",
                peak_reads.len(),
                real as u8,
                cap_signal as u8,
                reads[backbone_idx].name,
                windows.join(",")
            ));
        }
        if real {
            validated.push(peak);
        }
    }
    validated
}

/// What start validation sees across a library: its labelled reads (start more than
/// `START_RADIUS` into the backbone, with an exact match there), how many of them carry an
/// untemplated G, the capped share and sharpness of every start candidate with `min_iso`
/// labelled reads, and the sharpness of every start candidate.
#[derive(Default)]
pub(crate) struct StartStats {
    pub labelled: u64,
    pub capped: u64,
    /// (capped share, sharpness) of the candidates with `min_iso` labelled reads.
    pub candidates: Vec<(f32, f32)>,
    pub sharpness: Vec<f32>,
}

impl StartStats {
    pub(crate) fn merge(&mut self, other: StartStats) {
        self.labelled += other.labelled;
        self.capped += other.capped;
        self.candidates.extend(other.candidates);
        self.sharpness.extend(other.sharpness);
    }

    /// The library's cap signal and thresholds, unless set on the command line. The cap signal
    /// is on when at least `CAP_LIBRARY` of the labelled reads are capped, and, with at least
    /// `MIN_LIBRARY_CANDIDATES` labelled candidates, the candidates above the cap threshold are
    /// sharper than those below it (one-sided Mann-Whitney, p < 0.001): a cap marks real
    /// starts, which form tight piles. In ONT PCR-cDNA, untemplated-looking Gs left by primer
    /// trimming pass the share but not this test. The cap threshold is
    /// the cut between the two modes of the candidates' capped shares, truncation peaks near
    /// none and real starts near all (Otsu); the sharpness threshold is the Otsu cut over the
    /// candidates' sharpness, diffuse truncation against tight starts. With fewer than
    /// `MIN_LIBRARY_CANDIDATES` candidates, `MIN_CAPPED` and `MIN_SHARPNESS`.
    pub(crate) fn library(&self, options: &IsoformOptions) -> StartLibrary {
        let threshold = options.cap_threshold.unwrap_or_else(|| {
            if self.candidates.len() >= MIN_LIBRARY_CANDIDATES {
                let shares: Vec<f32> = self.candidates.iter().map(|&(share, _)| share).collect();
                otsu_threshold(&shares)
            } else {
                MIN_CAPPED
            }
        });
        let cap_signal = match options.cap_signal {
            CapSignal::On => true,
            CapSignal::Off => false,
            CapSignal::Auto => {
                self.labelled >= MIN_CAP_READS as u64
                    && self.capped as f64 >= CAP_LIBRARY * self.labelled as f64
                    && (self.candidates.len() < MIN_LIBRARY_CANDIDATES
                        || self.capped_starts_are_sharper(threshold))
            }
        };
        let sharpness = options.sharp_threshold.unwrap_or_else(|| {
            if self.sharpness.len() >= MIN_LIBRARY_CANDIDATES {
                otsu_threshold(&self.sharpness)
            } else {
                MIN_SHARPNESS
            }
        });
        StartLibrary { cap_signal, threshold, sharpness }
    }

    /// Whether the candidates with at least `threshold` capped reads are sharper than the
    /// others: one-sided Mann-Whitney U test, normal approximation, z > 3.09 (p < 0.001).
    pub(crate) fn capped_starts_are_sharper(&self, threshold: f64) -> bool {
        let (z, _) = self.cap_sharpness_test(threshold);
        z > 3.09
    }

    /// The Mann-Whitney z and AUC of sharpness, capped candidates against uncapped ones.
    pub(crate) fn cap_sharpness_test(&self, threshold: f64) -> (f64, f64) {
        let mut ranked: Vec<(f32, bool)> = self
            .candidates
            .iter()
            .map(|&(share, sharp)| (sharp, share as f64 >= threshold))
            .collect();
        ranked.sort_by(|a, b| a.0.total_cmp(&b.0));
        let (n1, n2) = ranked.iter().fold((0f64, 0f64), |(a, b), &(_, c)| {
            if c {
                (a + 1.0, b)
            } else {
                (a, b + 1.0)
            }
        });
        if n1 == 0.0 || n2 == 0.0 {
            return (0.0, 0.5);
        }
        // Sum of the capped candidates' ranks, ties sharing their mean rank.
        let mut rank_sum = 0f64;
        let mut i = 0;
        while i < ranked.len() {
            let mut j = i;
            while j < ranked.len() && ranked[j].0 == ranked[i].0 {
                j += 1;
            }
            let mean_rank = (i + 1 + j) as f64 / 2.0;
            rank_sum += mean_rank * ranked[i..j].iter().filter(|&&(_, c)| c).count() as f64;
            i = j;
        }
        let u = rank_sum - n1 * (n1 + 1.0) / 2.0;
        let z = (u - n1 * n2 / 2.0) / (n1 * n2 * (n1 + n2 + 1.0) / 12.0).sqrt();
        (z, u / (n1 * n2))
    }
}

/// Add what start validation sees in `group`, sitting on `reads[backbone_idx]`, to `acc`.
pub(super) fn add_start_stats(
    group: &[PlacedRead],
    reads: &[Read],
    backbone_idx: usize,
    options: &IsoformOptions,
    acc: &mut StartStats,
) {
    let width = options.peak_width.max(1) as i32;
    let backbone_seq = &reads[backbone_idx].seq;
    let capped: Vec<Option<bool>> = group
        .iter()
        .map(|&(read_idx, start, _)| {
            (start > START_RADIUS)
                .then(|| untemplated_lead(&reads[read_idx].seq, backbone_seq, start as usize))
                .flatten()
                .map(|(n, first)| n > 0 && first == b'G')
        })
        .collect();
    acc.labelled += capped.iter().flatten().count() as u64;
    acc.capped += capped.iter().flatten().filter(|&&c| c).count() as u64;
    let peaks = find_peaks(group.iter().map(|&(_, start, _)| start), width, options.min_iso);
    for peak in peaks.into_iter().filter(|&peak| peak > SITE_WIDTH) {
        let near: Vec<usize> = (0..group.len())
            .filter(|&group_pos| (group[group_pos].1 - peak).abs() <= START_RADIUS)
            .collect();
        let (n, g) = near
            .iter()
            .filter_map(|&group_pos| capped[group_pos])
            .fold((0usize, 0usize), |(n, g), c| (n + 1, g + c as usize));
        let around =
            group.iter().filter(|&&(_, start, _)| (start - peak).abs() <= SHARP_WINDOW).count();
        let sharpness = near.len() as f32 / around as f32;
        if n >= options.min_iso {
            acc.candidates.push((g as f32 / n as f32, sharpness));
        }
        acc.sharpness.push(sharpness);
    }
}

/// Otsu's threshold for values in [0, 1]: over 100 bins, the cut that maximises the variance
/// between the values below and above it.
fn otsu_threshold(values: &[f32]) -> f64 {
    const BINS: usize = 100;
    let mut hist = [0u64; BINS];
    for &v in values {
        hist[((v as f64 * BINS as f64) as usize).min(BINS - 1)] += 1;
    }
    let mid = |bin: usize| (bin as f64 + 0.5) / BINS as f64;
    let n: u64 = hist.iter().sum();
    let total: f64 = hist.iter().enumerate().map(|(bin, &c)| mid(bin) * c as f64).sum();
    let (mut below, mut below_sum) = (0u64, 0f64);
    let (mut best, mut cut) = (-1f64, MIN_CAPPED);
    for t in 1..BINS {
        below += hist[t - 1];
        below_sum += mid(t - 1) * hist[t - 1] as f64;
        let above = n - below;
        if below == 0 || above == 0 {
            continue;
        }
        let gap = below_sum / below as f64 - (total - below_sum) / above as f64;
        let between = below as f64 * above as f64 * gap * gap;
        if between > best {
            best = between;
            cut = t as f64 / BINS as f64;
        }
    }
    cut
}

/// How `read_seq`, whose first matched base sits at backbone position `start`, begins: the
/// number of leading bases (0-3) before it matches the backbone exactly for 12 bases there
/// (within 1 bp), and its first base. An untemplated G is 1-3 leading bases, the first a G,
/// that the backbone doesn't have there. None if the read's start doesn't match the backbone
/// exactly nearby (sequencing errors), so the read tells nothing. A G that the transcript
/// itself has just before the start is indistinguishable and counts as templated.
fn untemplated_lead(read_seq: &[u8], backbone_seq: &[u8], start: usize) -> Option<(usize, u8)> {
    const L: usize = 12;
    for lead in 0..=3usize {
        for b in [start, start.saturating_sub(1), start + 1] {
            if b + L <= backbone_seq.len()
                && lead + L <= read_seq.len()
                && read_seq[lead..lead + L] == backbone_seq[b..b + L]
            {
                return Some((lead, read_seq[0]));
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

    #[test]
    fn a_start_pile_within_a_site_width_of_the_backbone_start_is_no_candidate() {
        // Three reads at the backbone's start and six starting together `p` bp downstream
        // (validated as a sharp peak): a start of its own only beyond a site width.
        let alt_starts = |p: i32| {
            let mut group: Vec<PlacedRead> = (0..3).map(|i| (i, 0, 2000)).collect();
            group.extend((3..9).map(|i| (i, p, 2000)));
            let reads: Vec<Read> =
                (0..9).map(|i| Read { name: i.to_string(), seq: Vec::new() }).collect();
            let options = IsoformOptions { start_split: StartSplit::Auto, ..IsoformOptions::default() };
            split_by_ends(&group, &reads, 0, &options).iter().filter(|i| i.1).count()
        };
        assert_eq!(alt_starts(40), 0);
        assert_eq!(alt_starts(60), 1);
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
        let cfg = IsoformOptions { start_split: StartSplit::On, ..IsoformOptions::default() };
        assert_eq!(sizes(&truncated_majority(), &cfg), vec![8]);
    }

    #[test]
    fn without_start_split_full_length_reads_join_their_isoform() {
        let options = IsoformOptions { start_split: StartSplit::Off, ..IsoformOptions::default() };
        assert_eq!(sizes(&truncated_majority(), &options), vec![11]);
    }

    #[test]
    fn otsu_cuts_between_two_modes() {
        let mut values: Vec<f32> = (0..300).map(|i| (i % 10) as f32 / 100.0).collect();
        values.extend((0..100).map(|i| 0.75 + (i % 20) as f32 / 100.0));
        let cut = otsu_threshold(&values);
        assert!(cut > 0.09 && cut <= 0.75, "{cut}");
    }

    #[test]
    fn the_library_decides_the_cap_signal_and_threshold() {
        // capped candidates (share 0.9) are sharp, uncapped ones (0.05) diffuse
        let mut candidates: Vec<(f32, f32)> = (0..150).map(|i| (0.05, (i % 30) as f32 / 100.0)).collect();
        candidates.extend((0..50).map(|i| (0.9, 0.5 + (i % 40) as f32 / 100.0)));
        let stats = StartStats { labelled: 1000, capped: 400, candidates, sharpness: Vec::new() };
        let auto = IsoformOptions::default();
        assert!(stats.library(&auto).cap_signal);
        // the same Gs on candidates no sharper than the rest: not a cap signal (ONT PCR-cDNA)
        let mut mixed: Vec<(f32, f32)> = (0..150).map(|i| (0.05, (i % 50) as f32 / 100.0)).collect();
        mixed.extend((0..50).map(|i| (0.9, (i % 50) as f32 / 100.0)));
        let primer = StartStats { labelled: 1000, capped: 400, candidates: mixed, sharpness: Vec::new() };
        assert!(!primer.library(&auto).cap_signal);
        let few = StartStats { labelled: 1000, capped: 100, ..StartStats::default() };
        assert!(!few.library(&auto).cap_signal);
        assert_eq!(few.library(&auto).threshold, MIN_CAPPED);
        assert_eq!(few.library(&auto).sharpness, MIN_SHARPNESS);
        let set = IsoformOptions {
            cap_signal: CapSignal::On,
            cap_threshold: Some(0.5),
            sharp_threshold: Some(0.3),
            ..auto
        };
        assert_eq!(
            few.library(&set),
            StartLibrary { cap_signal: true, threshold: 0.5, sharpness: 0.3 }
        );
        // without a cap signal, the sharpness threshold comes from the candidates
        let mut sharp: Vec<f32> = (0..300).map(|i| (i % 20) as f32 / 100.0).collect();
        sharp.extend((0..100).map(|i| 0.8 + (i % 20) as f32 / 100.0));
        let uncapped = StartStats { labelled: 1000, capped: 50, sharpness: sharp, ..StartStats::default() };
        let cut = uncapped.library(&auto).sharpness;
        assert!((0.19..=0.8).contains(&cut), "{cut}");
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
