//! Step 4: merge isoforms that turn out to be the same transcript.
//!
//! Detection can split one transcript into several isoforms: reads whose errors break
//! their anchors fail step 1 and gather in small satellite groups, and a read rejected by
//! one backbone can found a duplicate group. Consensuses are error-corrected, so they can be
//! compared more strictly than reads were.
//!
//! Isoforms are visited largest first, and each merges into the first larger one kept so far
//! that either
//!   * has a near-identical consensus (`collapse_ident`), whatever the read counts, or
//!   * contains it by the step-1 test with `max_gap = collapse_gap` and `min_cov >= 0.95`,
//!     if it has at most `collapse_ratio` of that isoform's reads. Comparable read counts
//!     suggest real isoforms with different ends rather than a satellite.
//!
//! A merged isoform's reads join the kept isoform, which keeps its own consensus, unless the
//! merged one is a more complete copy of it — reaching further 5' than a 5'-truncated
//! majority, or on to the real 3' end past reads that stop early (internal priming): then the
//! kept isoform's consensus is rebuilt from all its reads, so the full-length sequence survives.
//! Containment treats an indel of more than `collapse_gap` bp as a difference, as the
//! near-identity test does.
//!
//! An isoform at a validated start is a transcript of its own, so a pair whose starts differ
//! stays apart when either is one. Two validated starts differ beyond `collapse_gap`; a
//! validated start and one that isn't (its backbone's start, the upstream edge of its start
//! site) differ only beyond a site width (`SITE_WIDTH`).
//!
//! Two switches change containment:
//!   * `--alt-starts keep` (the default): an isoform at a validated start never merges into
//!     one it lies inside with a different start; with `merge`, only the larger of a pair is
//!     protected;
//!   * `--end-merge polya`: a pair differing only at the 3' end is decided by internal
//!     priming instead of read counts. HiFi FLNC reads end at a polyA tail, so an earlier 3'
//!     end is either a real polyA site or the oligo-dT primer annealed to an A-rich stretch;
//!     the longer isoform holds the sequence just past that end, so an A-rich window there
//!     merges the pair and anything else keeps both.

use std::cmp::Reverse;
use std::collections::HashMap;

use super::anchors::{minimizer_map, MinimizerMap};
use super::ends::SITE_WIDTH;
use super::fit::{place_on_backbone, strip_homopolymer_ends};
use super::{align, consensus, AltStarts, EndMerge, Isoform, IsoformOptions, Read};

/// Merge duplicate isoforms among `isoforms`, whose members index `reads` (with minimizer maps
/// `read_maps`, used to rebuild a consensus).
pub(super) fn merge_duplicates(
    mut isoforms: Vec<Isoform>,
    reads: &[Read],
    read_maps: &[MinimizerMap],
    options: &IsoformOptions,
) -> Vec<Isoform> {
    if isoforms.len() < 2 {
        return isoforms;
    }
    isoforms.sort_by_key(|isoform| Reverse(isoform.members.len()));
    let consensus_maps: Vec<MinimizerMap> = isoforms
        .iter()
        .map(|isoform| minimizer_map(&isoform.consensus, options.minimizer_k, options.minimizer_w))
        .collect();
    let strict = IsoformOptions {
        min_cov: options.min_cov.max(0.95),
        max_gap: options.collapse_gap,
        max_indel_run: options.collapse_gap + 1,
        ..options.clone()
    };

    let n = isoforms.len();
    let mut merged_into: Vec<Option<usize>> = vec![None; n]; // merged isoform -> kept isoform
    let mut rebuild = vec![false; n]; // kept isoform absorbed a more complete copy of itself
    let mut kept: Vec<usize> = Vec::new();
    for i in 0..n {
        let mut into = None;
        for &j in &kept {
            match decide_merge(
                &isoforms[i],
                &isoforms[j],
                &consensus_maps[i],
                &consensus_maps[j],
                options,
                &strict,
            ) {
                MergeDecision::KeepSeparate => continue,
                MergeDecision::Merge => {}
                MergeDecision::MergeAndRebuild => rebuild[j] = true,
            }
            into = Some(j);
            break;
        }
        match into {
            Some(j) => merged_into[i] = Some(j),
            None => kept.push(i),
        }
    }

    // Kept isoforms in visiting order, each absorbing the reads merged into it.
    let mut pos: HashMap<usize, usize> = HashMap::new();
    let mut kept_isoforms: Vec<Isoform> = Vec::with_capacity(kept.len());
    for (p, &kept_idx) in kept.iter().enumerate() {
        pos.insert(kept_idx, p);
        kept_isoforms.push(Isoform {
            members: std::mem::take(&mut isoforms[kept_idx].members),
            consensus: std::mem::take(&mut isoforms[kept_idx].consensus),
            alt_start: isoforms[kept_idx].alt_start,
        });
    }
    for i in 0..n {
        if let Some(j) = merged_into[i] {
            let extra = std::mem::take(&mut isoforms[i].members);
            kept_isoforms[pos[&j]].members.extend(extra);
        }
    }
    for (p, &kept_idx) in kept.iter().enumerate() {
        if rebuild[kept_idx] {
            let isoform = &mut kept_isoforms[p];
            isoform.consensus =
                consensus::build_consensus(reads, &isoform.members, &isoform.members, read_maps);
        }
    }
    kept_isoforms
}

/// How `small` relates to a larger isoform it may merge into.
enum MergeDecision {
    KeepSeparate,
    Merge,
    /// Merge, and `small` is a more complete copy of `big`: rebuild the consensus.
    MergeAndRebuild,
}

/// Whether `small` (with no more reads than `big`) should merge into `big`. `strict` is
/// `options` tightened for comparing consensuses.
fn decide_merge(
    small: &Isoform,
    big: &Isoform,
    small_map: &MinimizerMap,
    big_map: &MinimizerMap,
    options: &IsoformOptions,
    strict: &IsoformOptions,
) -> MergeDecision {
    if consensuses_near_identical(small, big, options) {
        return MergeDecision::Merge;
    }
    // With `EndMerge::Polya`, a pair differing only at the 3' end is decided by internal
    // priming whatever the read counts, so containment is tested even past the ratio.
    let polya = options.end_merge == EndMerge::Polya;
    let within_ratio =
        small.members.len() as f64 <= options.collapse_ratio * big.members.len() as f64;
    if !within_ratio && !polya {
        return MergeDecision::KeepSeparate;
    }
    // Containment, testing the shorter consensus against the longer. An end reaching more
    // than `collapse_gap` bp past the other's is a real extension, not a consensus wobble.
    // An isoform at a validated downstream start is a transcript of its own, so a pair whose
    // starts differ stays apart when either is one: two validated starts beyond
    // `collapse_gap`, a validated start and its backbone's start beyond `SITE_WIDTH`.
    let end_slack = options.collapse_gap as i32;
    let start_slack = match (small.alt_start, big.alt_start) {
        (false, false) => None,
        (true, true) => Some(end_slack),
        _ => Some(SITE_WIDTH),
    };
    let starts_differ = |start: i32| start_slack.is_some_and(|slack| start.abs() > slack);
    if small.consensus.len() <= big.consensus.len() {
        // `small` placed on `big`; it may overhang `big`'s 5' end (a negative start).
        match place_on_backbone(small_map, big_map, &small.consensus, &big.consensus, strict) {
            None => MergeDecision::KeepSeparate,
            Some(_) if !within_ratio && !polya => MergeDecision::KeepSeparate,
            Some((start, end)) if start < -end_slack => {
                // Same 3' end: `small` is `big` with its 5' end intact. Otherwise the two
                // differ at both ends and stay apart.
                if within_ratio
                    && !starts_differ(start)
                    && end >= big.consensus.len() as i32 - options.boundary_tol as i32
                {
                    MergeDecision::MergeAndRebuild
                } else {
                    MergeDecision::KeepSeparate
                }
            }
            Some((start, _))
                if options.alt_starts == AltStarts::Keep && small.alt_start && starts_differ(start) =>
            {
                MergeDecision::KeepSeparate
            }
            Some((_, end)) if polya && end < big.consensus.len() as i32 - end_slack => {
                // `small` ends early: internal priming if `big` is A-rich just past its end.
                if internally_primed(&big.consensus, end, options) {
                    MergeDecision::Merge
                } else {
                    MergeDecision::KeepSeparate
                }
            }
            Some(_) if within_ratio => MergeDecision::Merge,
            Some(_) => MergeDecision::KeepSeparate,
        }
    } else {
        // `big` placed on `small`: `small` extends it at the 5' end, the 3' end, or both.
        match place_on_backbone(big_map, small_map, &big.consensus, &small.consensus, strict) {
            None => MergeDecision::KeepSeparate,
            Some((start, _)) if starts_differ(start) => MergeDecision::KeepSeparate,
            Some((_, end)) if polya && end < small.consensus.len() as i32 - end_slack => {
                // `big` ends early: internal priming if `small` is A-rich just past its end.
                if internally_primed(&small.consensus, end, options) {
                    MergeDecision::MergeAndRebuild
                } else {
                    MergeDecision::KeepSeparate
                }
            }
            Some(_) if !within_ratio => MergeDecision::KeepSeparate,
            Some((start, end))
                if start > end_slack || end < small.consensus.len() as i32 - end_slack =>
            {
                MergeDecision::MergeAndRebuild
            }
            Some(_) => MergeDecision::Merge,
        }
    }
}

/// Whether the 3' end at `end` of a shorter isoform, placed on `longer`, looks like internal
/// priming: the oligo-dT primer annealed to an A-rich stretch of the transcript, which
/// `longer` holds just past `end`.
fn internally_primed(longer: &[u8], end: i32, options: &IsoformOptions) -> bool {
    let from = (end.max(0) as usize).min(longer.len());
    let window = &longer[from..(from + options.priming_window).min(longer.len())];
    if window.is_empty() {
        return false;
    }
    let a = window.iter().filter(|&&b| b == b'A' || b == b'a').count();
    a as f64 >= options.priming_a_frac * window.len() as f64
}

/// Near-identical consensuses, ignoring homopolymer ends: identity at least `collapse_ident`,
/// with the longest indel and the length difference both within `collapse_gap`.
fn consensuses_near_identical(a: &Isoform, b: &Isoform, options: &IsoformOptions) -> bool {
    let a_seq = strip_homopolymer_ends(&a.consensus);
    let b_seq = strip_homopolymer_ends(&b.consensus);
    let alignment = align::align_global(a_seq, b_seq);
    alignment.ident() >= options.collapse_ident
        && alignment.longest_indel <= options.collapse_gap as usize
        && (a_seq.len() as i64 - b_seq.len() as i64).unsigned_abs() as u32 <= options.collapse_gap
}

#[cfg(test)]
mod tests {
    use super::*;

    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn a_5_prime_complete_minority_gives_the_merged_isoform_its_full_sequence() {
        // Ten reads missing the first 500 bp form the larger isoform; three full-length reads
        // form a smaller one. They merge, and the kept isoform must not stay truncated.
        let truth = random_seq(1500, 5);
        let mut seqs: Vec<Vec<u8>> = (0..10).map(|_| truth[500..].to_vec()).collect();
        seqs.extend((0..3).map(|_| truth.clone()));
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(read_idx, seq)| Read { name: read_idx.to_string(), seq: seq.clone() })
            .collect();
        let options = IsoformOptions::default();
        let read_maps: Vec<MinimizerMap> = reads
            .iter()
            .map(|r| minimizer_map(&r.seq, options.minimizer_k, options.minimizer_w))
            .collect();
        let isoforms = vec![
            Isoform {
                members: (0..10).collect(),
                consensus: truth[500..].to_vec(),
                alt_start: false,
            },
            Isoform {
                members: (10..13).collect(),
                consensus: truth.clone(),
                alt_start: false,
            },
        ];
        let out = merge_duplicates(isoforms, &reads, &read_maps, &options);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].members.len(), 13);
        assert_eq!(out[0].consensus, truth);
    }

    #[test]
    fn an_isoform_at_a_validated_start_is_not_folded_into_the_full_length_one() {
        let truth = random_seq(1500, 5);
        let mut seqs: Vec<Vec<u8>> = (0..10).map(|_| truth[500..].to_vec()).collect();
        seqs.extend((0..3).map(|_| truth.clone()));
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(read_idx, seq)| Read { name: read_idx.to_string(), seq: seq.clone() })
            .collect();
        let options = IsoformOptions::default();
        let read_maps: Vec<MinimizerMap> = reads
            .iter()
            .map(|r| minimizer_map(&r.seq, options.minimizer_k, options.minimizer_w))
            .collect();
        let isoforms = vec![
            Isoform {
                members: (0..10).collect(),
                consensus: truth[500..].to_vec(),
                alt_start: true,
            },
            Isoform {
                members: (10..13).collect(),
                consensus: truth.clone(),
                alt_start: false,
            },
        ];
        assert_eq!(merge_duplicates(isoforms, &reads, &read_maps, &options).len(), 2);
    }

    #[test]
    fn reads_reaching_the_real_3_prime_end_give_the_merged_isoform_its_full_sequence() {
        // Ten reads stop 600 bp early (internal priming) and form the larger isoform; three
        // reach the real 3' end. The merged isoform must keep the full-length sequence.
        let truth = random_seq(1500, 9);
        let mut seqs: Vec<Vec<u8>> = (0..10).map(|_| truth[..900].to_vec()).collect();
        seqs.extend((0..3).map(|_| truth.clone()));
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(read_idx, seq)| Read { name: read_idx.to_string(), seq: seq.clone() })
            .collect();
        let options = IsoformOptions::default();
        let read_maps: Vec<MinimizerMap> = reads
            .iter()
            .map(|r| minimizer_map(&r.seq, options.minimizer_k, options.minimizer_w))
            .collect();
        let isoforms = vec![
            Isoform {
                members: (0..10).collect(),
                consensus: truth[..900].to_vec(),
                alt_start: false,
            },
            Isoform {
                members: (10..13).collect(),
                consensus: truth.clone(),
                alt_start: false,
            },
        ];
        let out = merge_duplicates(isoforms, &reads, &read_maps, &options);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].consensus, truth);
    }

    /// Isoforms of `counts[i]` identical reads with consensus `seqs[i]`, run through merging;
    /// returns (reads, consensus length) per final isoform, largest first.
    fn merge_case(seqs: &[Vec<u8>], counts: &[usize], alt: &[bool], options: &IsoformOptions)
        -> Vec<(usize, usize)> {
        let mut all: Vec<Vec<u8>> = Vec::new();
        let mut isoforms = Vec::new();
        for ((seq, &n), &a) in seqs.iter().zip(counts).zip(alt) {
            let first = all.len();
            all.extend((0..n).map(|_| seq.clone()));
            isoforms.push(Isoform { members: (first..first + n).collect(), consensus: seq.clone(), alt_start: a });
        }
        let reads: Vec<Read> = all.iter().enumerate()
            .map(|(i, seq)| Read { name: i.to_string(), seq: seq.clone() }).collect();
        let read_maps: Vec<MinimizerMap> = reads.iter()
            .map(|r| minimizer_map(&r.seq, options.minimizer_k, options.minimizer_w)).collect();
        merge_duplicates(isoforms, &reads, &read_maps, options)
            .iter().map(|iso| (iso.members.len(), iso.consensus.len())).collect()
    }

    #[test]
    fn a_smaller_isoform_at_a_validated_start_is_kept_with_alt_starts_keep() {
        // 12 full-length reads and 4 starting at a validated start 500 bp downstream: with
        // --alt-starts merge the smaller one merges by containment; with keep it stays.
        let truth = random_seq(1500, 5);
        let seqs = [truth.clone(), truth[500..].to_vec()];
        let merge = IsoformOptions { alt_starts: AltStarts::Merge, ..IsoformOptions::default() };
        assert_eq!(merge_case(&seqs, &[12, 4], &[false, true], &merge), vec![(16, 1500)]);
        let keep = IsoformOptions { alt_starts: AltStarts::Keep, ..IsoformOptions::default() };
        assert_eq!(merge_case(&seqs, &[12, 4], &[false, true], &keep), vec![(12, 1500), (4, 1000)]);
        // a truncated copy without a validated start still merges
        assert_eq!(merge_case(&seqs, &[12, 4], &[false, false], &keep), vec![(16, 1500)]);
    }

    #[test]
    fn a_validated_start_within_a_site_width_of_one_that_is_not_validated_is_the_same_start() {
        // A full-length isoform and one at a validated start `d` bp downstream: within a site
        // width they are one start and merge; further apart, or both validated, they stay.
        let truth = random_seq(1500, 5);
        let at = |d: usize| [truth.clone(), truth[d..].to_vec()];
        let keep = IsoformOptions::default();
        assert_eq!(merge_case(&at(40), &[12, 4], &[false, true], &keep), vec![(16, 1500)]);
        assert_eq!(merge_case(&at(80), &[12, 4], &[false, true], &keep).len(), 2);
        assert_eq!(merge_case(&at(40), &[12, 4], &[true, true], &keep).len(), 2);
        // the full-length one the smaller: it extends the validated one, which is rebuilt
        let longer_minor = [truth[40..].to_vec(), truth.clone()];
        assert_eq!(merge_case(&longer_minor, &[12, 4], &[true, false], &keep), vec![(16, 1500)]);
    }

    #[test]
    fn end_merge_polya_keeps_a_minor_real_3_prime_end() {
        // A proximal end at 1500 that the distal isoform continues past with ordinary sequence.
        let truth = random_seq(2000, 9);
        let seqs = [truth.clone(), truth[..1500].to_vec()];
        let ratio = IsoformOptions::default();
        assert_eq!(merge_case(&seqs, &[10, 4], &[false, false], &ratio), vec![(14, 2000)]);
        let polya = IsoformOptions { end_merge: EndMerge::Polya, ..IsoformOptions::default() };
        assert_eq!(merge_case(&seqs, &[10, 4], &[false, false], &polya), vec![(10, 2000), (4, 1500)]);
        // the same when the distal end is the minor one
        assert_eq!(merge_case(&seqs, &[4, 10], &[false, false], &polya), vec![(10, 1500), (4, 2000)]);
    }

    #[test]
    fn end_merge_polya_merges_an_internally_primed_end_whatever_the_read_counts() {
        // The distal isoform is A-rich just past 1500: reads ending there are internal priming
        // and merge even with comparable read counts, keeping the distal sequence.
        let mut truth = random_seq(2000, 9);
        for (i, b) in b"AAAAGAAAAAAAACAAAAAA".iter().enumerate() {
            truth[1500 + i] = *b;
        }
        let seqs = [truth.clone(), truth[..1500].to_vec()];
        let polya = IsoformOptions { end_merge: EndMerge::Polya, ..IsoformOptions::default() };
        assert_eq!(merge_case(&seqs, &[10, 8], &[false, false], &polya), vec![(18, 2000)]);
        assert_eq!(merge_case(&seqs, &[8, 10], &[false, false], &polya), vec![(18, 2000)]);
        let ratio = IsoformOptions::default();
        assert_eq!(merge_case(&seqs, &[10, 8], &[false, false], &ratio).len(), 2);
    }
}
