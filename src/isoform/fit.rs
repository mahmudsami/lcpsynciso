//! Step 1: does a read have the same structure as its group's backbone?
//!
//! Minimizer anchors only show where the two sequences agree. The regions between and
//! around them decide the answer, compared by alignment ([`align`]):
//!
//!   * interior gaps, between consecutive anchors: a gap wider than `max_gap` is aligned and
//!     must be the same sequence (`min_gap_ident`) without one long indel (`max_indel_run`);
//!   * flanks, before the first anchor and after the last: each is extended outward from its
//!     anchor. A flank left unmatched while backbone sequence remains means the sequences
//!     diverge; a flank unmatched only because the backbone ran out is an overhang.

use std::time::Instant;

use super::anchors::{self, Anchor, MinimizerMap};
use super::stats::{self, RejectReason};
use super::{align, IsoformOptions};

/// Test whether `read_seq` has the same structure as `backbone_seq`. If it does, return the
/// read's start and end in backbone coordinates; a read that overhangs the backbone gets a
/// negative start or an end beyond the backbone's length.
///
/// `read_map` and `backbone_map` are the minimizer maps of the two sequences.
pub(super) fn place_on_backbone(
    read_map: &MinimizerMap,
    backbone_map: &MinimizerMap,
    read_seq: &[u8],
    backbone_seq: &[u8],
    options: &IsoformOptions,
) -> Option<(i32, i32)> {
    stats::inc(&stats::N_PLACEMENT_CALLS);

    // ── Anchors ──
    // Every shared minimizer k-mer pairs each of its read positions with each of its backbone
    // positions; a k-mer unique to both sequences (the common case) contributes exactly one.
    let t_anchor = Instant::now();
    let mut candidate_anchors = anchors::shared_anchors(read_map, backbone_map);
    if candidate_anchors.len() < options.min_anchors {
        stats::add_elapsed(&stats::T_ANCHOR, t_anchor);
        return reject(RejectReason::FewAnchors);
    }
    let chain = anchors::chain_anchors(&mut candidate_anchors);
    stats::add_elapsed(&stats::T_ANCHOR, t_anchor);
    if chain.len() < options.min_anchors {
        return reject(RejectReason::ShortChain);
    }

    // ── Interior gaps ──
    let t_gap = Instant::now();
    let gaps_match = interior_gaps_match(&chain, read_seq, backbone_seq, options);
    stats::add_elapsed(&stats::T_GAP, t_gap);
    if !gaps_match {
        return reject(RejectReason::InteriorGap);
    }

    // ── Flanks, each extended outward from its terminal anchor ──
    let t_flank = Instant::now();
    let (first_read_pos, first_backbone_pos) = chain[0];
    let (last_read_pos, last_backbone_pos) = *chain.last().unwrap();

    // 5' flank, aligned reversed so both sequences start at the anchor.
    let head_read = &read_seq[..first_read_pos as usize];
    let head_backbone = &backbone_seq[..first_backbone_pos as usize];
    stats::record_flank_size(head_read.len().min(head_backbone.len()));
    let head_fit = align::extend_flank(head_read, head_backbone, true);
    let start = if head_fit.backbone_exhausted {
        // The read starts before the backbone, unless that overhang is a homopolymer.
        let over = head_read.len().saturating_sub(head_fit.backbone_consumed);
        let unmatched = head_fit.unmatched.min(head_read.len());
        if unmatched > 0 && is_homopolymer_tail(&head_read[..unmatched]) {
            0
        } else {
            -(over as i32)
        }
    } else {
        first_backbone_pos as i32 - head_fit.backbone_consumed as i32
    };

    // 3' flank.
    let tail_backbone_start =
        (last_backbone_pos as usize + options.minimizer_k).min(backbone_seq.len());
    let tail_read = &read_seq[(last_read_pos as usize + options.minimizer_k).min(read_seq.len())..];
    let tail_backbone = &backbone_seq[tail_backbone_start..];
    stats::record_flank_size(tail_read.len().min(tail_backbone.len()));
    let tail_fit = align::extend_flank(tail_read, tail_backbone, false);
    let end = if tail_fit.backbone_exhausted {
        // The read ends after the backbone, unless that overhang is a homopolymer.
        let over = tail_read.len().saturating_sub(tail_fit.backbone_consumed);
        let unmatched = tail_fit.unmatched.min(tail_read.len());
        if unmatched > 0 && is_homopolymer_tail(&tail_read[tail_read.len() - unmatched..]) {
            (tail_backbone_start + tail_backbone.len()) as i32
        } else {
            (tail_backbone_start + tail_backbone.len()) as i32 + over as i32
        }
    } else {
        (tail_backbone_start + tail_fit.backbone_consumed) as i32
    };
    stats::add_elapsed(&stats::T_FLANK, t_flank);

    // ── Decide from the unmatched flanks ──
    // A homopolymer overhang never counts as unmatched.
    let head_is_homopolymer = head_fit.unmatched > 0
        && is_homopolymer_tail(&head_read[..head_fit.unmatched.min(head_read.len())]);
    let tail_is_homopolymer = tail_fit.unmatched > 0
        && is_homopolymer_tail(
            &tail_read[tail_read.len() - tail_fit.unmatched.min(tail_read.len())..],
        );
    let head_long = head_fit.unmatched as u32 > options.max_flank && !head_is_homopolymer;
    let tail_long = tail_fit.unmatched as u32 > options.max_flank && !tail_is_homopolymer;
    // Unmatched while backbone sequence remained: the sequences diverge at that end.
    if (head_long && !head_fit.backbone_exhausted) || (tail_long && !tail_fit.backbone_exhausted) {
        return reject(RejectReason::FlankDiverged);
    }
    // Overhangs at both ends: the read shares only its middle with the backbone.
    if head_long && tail_long {
        return reject(RejectReason::BothFlanks);
    }

    // ── Coverage ──
    let charged_head = if head_is_homopolymer { 0 } else { head_fit.unmatched };
    let charged_tail = if tail_is_homopolymer { 0 } else { tail_fit.unmatched };
    let covered = read_seq.len().saturating_sub(charged_head + charged_tail);
    if (covered as f64) / (read_seq.len().max(1) as f64) < options.min_cov {
        return reject(RejectReason::LowCoverage);
    }
    Some((start, end))
}

fn reject(reason: RejectReason) -> Option<(i32, i32)> {
    stats::count_rejection(reason);
    None
}

/// Check each gap between consecutive chain anchors. A gap of at most `max_gap` bases on
/// both sides is accepted as sequencing error; a wider one is aligned, and must be the same
/// sequence without one long indel.
fn interior_gaps_match(
    chain: &[Anchor],
    read_seq: &[u8],
    backbone_seq: &[u8],
    options: &IsoformOptions,
) -> bool {
    for pair in chain.windows(2) {
        let (prev_read_pos, prev_backbone_pos) = pair[0];
        let (next_read_pos, next_backbone_pos) = pair[1];
        // A gap runs from the end of one anchor's k-mer to the start of the next anchor.
        let read_gap_start = (prev_read_pos as usize + options.minimizer_k).min(read_seq.len());
        let read_gap_end = (next_read_pos as usize).min(read_seq.len());
        let backbone_gap_start =
            (prev_backbone_pos as usize + options.minimizer_k).min(backbone_seq.len());
        let backbone_gap_end = (next_backbone_pos as usize).min(backbone_seq.len());
        let gap_len = read_gap_end
            .saturating_sub(read_gap_start)
            .max(backbone_gap_end.saturating_sub(backbone_gap_start));
        if gap_len <= options.max_gap as usize {
            continue;
        }
        stats::inc(&stats::N_GAPS_ALIGNED);
        stats::record_gap_size(gap_len);
        stats::inc(&stats::N_TRACE);
        let gap_stats = align::align_global(
            &read_seq[read_gap_start..read_gap_end.max(read_gap_start)],
            &backbone_seq[backbone_gap_start..backbone_gap_end.max(backbone_gap_start)],
        );
        if gap_stats.ident() < options.min_gap_ident {
            return false; // different sequence
        }
        if gap_stats.longest_indel >= options.max_indel_run as usize {
            return false; // extra sequence inserted in one piece
        }
        stats::inc(&stats::N_GAPS_ACCEPTED);
    }
    true
}

// ── Homopolymer tails ────────────────────────────────────────────────────────

/// Shortest stretch (bp) that can count as a homopolymer tail.
const MIN_HOMOPOLYMER: usize = 10;
/// Fraction of a tail that must be a single base.
const HOMOPOLYMER_PURITY: f64 = 0.8;

/// Whether an unmatched read end is a homopolymer tail, such as polyA, rather than
/// transcript sequence.
fn is_homopolymer_tail(seq: &[u8]) -> bool {
    if seq.len() < MIN_HOMOPOLYMER {
        return false;
    }
    let mut counts = [0usize; 4];
    for &base in seq {
        match base {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {}
        }
    }
    *counts.iter().max().unwrap() as f64 / seq.len() as f64 >= HOMOPOLYMER_PURITY
}

/// `seq` with a run of more than `MIN_HOMOPOLYMER` identical bases at either end cut back to
/// one base, so two consensuses compare on transcript sequence rather than tail length.
pub(super) fn strip_homopolymer_ends(seq: &[u8]) -> &[u8] {
    let mut start = 0usize;
    let mut end = seq.len();
    while end > start + 1 && seq[end - 1].eq_ignore_ascii_case(&seq[end - 2]) {
        end -= 1;
    }
    if seq.len() - end < MIN_HOMOPOLYMER {
        end = seq.len();
    }
    while start + 1 < end && seq[start].eq_ignore_ascii_case(&seq[start + 1]) {
        start += 1;
    }
    if start < MIN_HOMOPOLYMER {
        start = 0;
    }
    &seq[start..end]
}
