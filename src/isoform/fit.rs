//! Step 1: does a read have the same structure as its group's backbone?
//!
//! Minimizer anchors only show where the two sequences agree. The regions between and
//! around them decide the answer, compared by alignment ([`ba`]):
//!
//!   * interior gaps, between consecutive anchors: a gap wider than `max_gap` is aligned and
//!     must be the same sequence (`min_gap_ident`) without one long indel (`max_indel_run`);
//!   * flanks, before the first anchor and after the last: each is extended outward from its
//!     anchor. A flank left unmatched while backbone sequence remains means the sequences
//!     diverge; a flank unmatched only because the backbone ran out is an overhang.

use std::collections::HashMap;
use std::time::Instant;

use super::stats::{self, Reject};
use super::{anchors, ba, Cfg};

/// Test whether `read` has the same structure as `backbone`. If it does, return the read's
/// start and end in backbone coordinates; a read that overhangs the backbone gets a negative
/// start or an end beyond the backbone's length.
///
/// `read_map` and `backbone_map` are the minimizer maps of the two sequences.
pub(super) fn fits_interval(
    read_map: &HashMap<u64, Vec<u32>>,
    backbone_map: &HashMap<u64, Vec<u32>>,
    read: &[u8],
    backbone: &[u8],
    cfg: &Cfg,
) -> Option<(i32, i32)> {
    stats::inc(&stats::N_CALLS);

    // ── Anchors ──
    // Every shared minimizer code pairs each of its read positions with each of its backbone
    // positions; a code unique to both sequences (the common case) contributes exactly one.
    let t_anchor = Instant::now();
    let mut anchors: Vec<(u32, u32)> = Vec::new();
    for (c, is) in read_map.iter() {
        if let Some(js) = backbone_map.get(c) {
            for &i in is {
                for &j in js {
                    anchors.push((i, j));
                }
            }
        }
    }
    if anchors.len() < cfg.min_anchors {
        stats::add_elapsed(&stats::T_ANCHOR, t_anchor);
        return reject(Reject::FewAnchors);
    }
    let chain = anchors::chain_anchors(&mut anchors);
    stats::add_elapsed(&stats::T_ANCHOR, t_anchor);
    if chain.len() < cfg.min_anchors {
        return reject(Reject::ShortChain);
    }

    // ── Interior gaps ──
    let t_gap = Instant::now();
    let gaps_match = interior_gaps_match(&chain, read, backbone, cfg);
    stats::add_elapsed(&stats::T_GAP, t_gap);
    if !gaps_match {
        return reject(Reject::InteriorGap);
    }

    // ── Flanks, each extended outward from its terminal anchor ──
    let t_flank = Instant::now();
    let (first_read, first_bb) = chain[0];
    let (last_read, last_bb) = *chain.last().unwrap();

    // 5' flank, aligned reversed so both sequences start at the anchor.
    let head_read = &read[..first_read as usize];
    let head_bb = &backbone[..first_bb as usize];
    stats::record_flank_size(head_read.len().min(head_bb.len()));
    let head = ba::flank(head_read, head_bb, true);
    let start = if head.reference_exhausted {
        // The read starts before the backbone, unless that overhang is a homopolymer.
        let over = head_read.len().saturating_sub(head.consumed);
        let unmatched = head.unmatched.min(head_read.len());
        if cfg.polya_clamp && unmatched > 0 && is_homopolymer_tail(&head_read[..unmatched]) {
            0
        } else {
            -(over as i32)
        }
    } else {
        first_bb as i32 - head.consumed as i32
    };

    // 3' flank.
    let tail_bb_start = (last_bb as usize + cfg.k).min(backbone.len());
    let tail_read = &read[(last_read as usize + cfg.k).min(read.len())..];
    let tail_bb = &backbone[tail_bb_start..];
    stats::record_flank_size(tail_read.len().min(tail_bb.len()));
    let tail = ba::flank(tail_read, tail_bb, false);
    let end = if tail.reference_exhausted {
        // The read ends after the backbone, unless that overhang is a homopolymer.
        let over = tail_read.len().saturating_sub(tail.consumed);
        let unmatched = tail.unmatched.min(tail_read.len());
        if cfg.polya_clamp
            && unmatched > 0
            && is_homopolymer_tail(&tail_read[tail_read.len() - unmatched..])
        {
            (tail_bb_start + tail_bb.len()) as i32
        } else {
            (tail_bb_start + tail_bb.len()) as i32 + over as i32
        }
    } else {
        (tail_bb_start + tail.consumed) as i32
    };
    stats::add_elapsed(&stats::T_FLANK, t_flank);

    // ── Decide from the unmatched flanks ──
    // With the clamp on, a homopolymer overhang never counts as unmatched.
    let head_poly = cfg.polya_clamp
        && head.unmatched > 0
        && is_homopolymer_tail(&head_read[..head.unmatched.min(head_read.len())]);
    let tail_poly = cfg.polya_clamp
        && tail.unmatched > 0
        && is_homopolymer_tail(&tail_read[tail_read.len() - tail.unmatched.min(tail_read.len())..]);
    let head_long = head.unmatched as u32 > cfg.max_flank && !head_poly;
    let tail_long = tail.unmatched as u32 > cfg.max_flank && !tail_poly;
    // Unmatched while backbone sequence remained: the sequences diverge at that end.
    if (head_long && !head.reference_exhausted) || (tail_long && !tail.reference_exhausted) {
        return reject(Reject::FlankDiverged);
    }
    // Overhangs at both ends: the read shares only its middle with the backbone.
    if head_long && tail_long {
        return reject(Reject::BothFlanks);
    }

    // ── Coverage ──
    let charged_head = if head_poly { 0 } else { head.unmatched };
    let charged_tail = if tail_poly { 0 } else { tail.unmatched };
    let covered = read.len().saturating_sub(charged_head + charged_tail);
    if (covered as f64) / (read.len().max(1) as f64) < cfg.min_cov {
        return reject(Reject::LowCoverage);
    }
    Some((start, end))
}

fn reject(why: Reject) -> Option<(i32, i32)> {
    stats::count_rejection(why);
    None
}

/// Check each gap between consecutive chain anchors. A gap of at most `max_gap` bases on
/// both sides is accepted as sequencing error; a wider one is aligned, and must be the same
/// sequence without one long indel.
fn interior_gaps_match(chain: &[(u32, u32)], read: &[u8], backbone: &[u8], cfg: &Cfg) -> bool {
    for pair in chain.windows(2) {
        // A gap runs from the end of one anchor's k-mer to the start of the next anchor.
        let r0 = (pair[0].0 as usize + cfg.k).min(read.len());
        let r1 = (pair[1].0 as usize).min(read.len());
        let b0 = (pair[0].1 as usize + cfg.k).min(backbone.len());
        let b1 = (pair[1].1 as usize).min(backbone.len());
        let gap_len = r1.saturating_sub(r0).max(b1.saturating_sub(b0));
        if gap_len <= cfg.max_gap as usize {
            continue;
        }
        stats::inc(&stats::N_GAPDP);
        stats::record_gap_size(gap_len);
        stats::inc(&stats::N_TRACE);
        let gap = ba::gap_runs(&read[r0..r1.max(r0)], &backbone[b0..b1.max(b0)]);
        if cfg.min_gap_ident > 0.0 && gap.ident() < cfg.min_gap_ident {
            return false; // different sequence
        }
        if cfg.max_indel_run > 0 && gap.max_indel >= cfg.max_indel_run as usize {
            return false; // extra sequence inserted in one piece
        }
        stats::inc(&stats::N_RESCUE);
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
fn is_homopolymer_tail(s: &[u8]) -> bool {
    if s.len() < MIN_HOMOPOLYMER {
        return false;
    }
    let mut counts = [0usize; 4];
    for &b in s {
        match b {
            b'A' | b'a' => counts[0] += 1,
            b'C' | b'c' => counts[1] += 1,
            b'G' | b'g' => counts[2] += 1,
            b'T' | b't' => counts[3] += 1,
            _ => {}
        }
    }
    *counts.iter().max().unwrap() as f64 / s.len() as f64 >= HOMOPOLYMER_PURITY
}

/// `s` with a run of more than `MIN_HOMOPOLYMER` identical bases at either end cut back to
/// one base, so two consensuses compare on transcript sequence rather than tail length.
pub(super) fn strip_homopolymer_ends(s: &[u8]) -> &[u8] {
    let mut a = 0usize;
    let mut b = s.len();
    while b > a + 1 && s[b - 1].eq_ignore_ascii_case(&s[b - 2]) {
        b -= 1;
    }
    if s.len() - b < MIN_HOMOPOLYMER {
        b = s.len();
    }
    while a + 1 < b && s[a].eq_ignore_ascii_case(&s[a + 1]) {
        a += 1;
    }
    if a < MIN_HOMOPOLYMER {
        a = 0;
    }
    &s[a..b]
}
