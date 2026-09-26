//! Isoform resolution: split the reads of one cluster into isoforms, each with a consensus
//! sequence. Used by `predict` after clustering, and by `find-isoforms` on an existing
//! cluster assignment.
//!
//! [`resolve_cluster`] runs five steps on a cluster:
//!
//!   1. Group by structure ([`fit`]). Reads are taken longest first: each round the longest
//!      unassigned read is the backbone, every other unassigned read that fits it joins its
//!      group, and the rest wait for a later round.
//!   2. Split a group where many reads share one indel at the same position ([`variant`];
//!      off by default). This runs before step 3, which would otherwise mix two structures
//!      that share their ends.
//!   3. Split each group by where its reads start and end ([`ends`]).
//!   4. Keep groups of at least `min_iso` reads and build each one's consensus ([`consensus`]).
//!   5. Merge isoforms whose consensuses show they are the same transcript ([`collapse`]).
//!
//! Which test tells two isoforms apart:
//!
//!   extra sequence in one piece (exon, retained intron)   interior-gap indel run   step 1
//!   same-length but different sequence (MXE)              interior-gap identity    step 1
//!   different sequence at one end                         unmatched flank          step 1
//!   a splice site shifted by a few bases                  recurrent indel          step 2
//!   different transcription start or end                  end peaks                step 3
//!     (with `--start-split off|auto` 3' ends only, plus validated starts for auto: 5'
//!     starts in cDNA mostly mark truncation)
//!   the same transcript found by two groups               consensus identity       step 5
//!
//! Differences confined to a flank shorter than `max_flank` are not detected. The defaults
//! were tuned on PacBio HiFi SIRV reads; ONT data has not been measured.

use std::cmp::Reverse;
use std::time::Instant;

mod anchors; // minimizer anchors and their collinear chain
mod ba; // block-aligner calls: flank extension, gap alignment, indel events
mod cli; // the `find-isoforms` subcommand
mod collapse; // step 5
mod consensus; // step 4
mod ends; // step 3
mod fit; // step 1
mod options; // Cfg, its defaults, flags and help text
mod output; // output files
pub(crate) mod stats; // diagnostic counters
mod variant; // step 2

use anchors::MinimizerMap;
pub use cli::{parse_args, run, usage, Config};
pub use options::{Cfg, StartSplit};
pub(crate) use options::{parse_flag, RESOLVE_HELP};
pub(crate) use output::{render, OutputFiles, Rendered};

pub struct Read {
    pub name: String,
    pub seq: Vec<u8>,
}

/// One resolved isoform.
pub(crate) struct Isoform {
    /// Member reads, as indices into the cluster's reads.
    pub members: Vec<usize>,
    /// Consensus sequence, or the longest member read when consensus is off.
    pub consensus: Vec<u8>,
    /// Starts at a downstream transcription start that `--start-split auto` validated, so
    /// collapse must not fold it into a longer isoform that differs only at the 5' end.
    pub alt_start: bool,
}

/// A read placed on its group's backbone: (read index, start, end), with start and end in
/// backbone coordinates.
type Member = (usize, i32, i32);

/// Resolve one cluster's reads into isoforms. Clusters are independent, so callers can
/// resolve many in parallel.
pub(crate) fn resolve_cluster(reads: &[Read], cfg: &Cfg) -> Vec<Isoform> {
    let t = Instant::now();
    let maps: Vec<MinimizerMap> =
        reads.iter().map(|r| anchors::minimizer_map(&r.seq, cfg.k, cfg.w)).collect();
    stats::add_elapsed(&stats::T_MAPS, t);

    let t = Instant::now();
    let (groups, parts) = partition(reads, &maps, cfg);
    stats::add_elapsed(&stats::T_SPLIT, t);

    let t = Instant::now();
    let isoforms: Vec<Isoform> = groups
        .into_iter()
        .map(|(members, alt_start, part)| {
            let consensus = if cfg.consensus {
                consensus::refine_consensus_with_maps(reads, &members, &parts[part], &maps)
            } else {
                let longest = *members.iter().max_by_key(|&&i| reads[i].seq.len()).unwrap();
                reads[longest].seq.clone()
            };
            Isoform { members, consensus, alt_start }
        })
        .collect();
    stats::add_elapsed(&stats::T_CONS, t);

    if cfg.collapse {
        collapse::collapse(isoforms, reads, &maps, cfg)
    } else {
        isoforms
    }
}

/// Steps 1–3 and the `min_iso` filter. Each isoform, largest first, as its member read
/// indices, whether it starts at a validated downstream start, and the index of its
/// structure group (after any variant split) in the second list, which holds each group's
/// reads: they share the isoform's structure, so they may vote in its consensus.
fn partition(
    reads: &[Read],
    maps: &[MinimizerMap],
    cfg: &Cfg,
) -> (Vec<(Vec<usize>, bool, usize)>, Vec<Vec<usize>>) {
    let mut isoforms: Vec<(Vec<usize>, bool, usize)> = Vec::new();
    let mut parts: Vec<Vec<usize>> = Vec::new();
    for group in &group_by_structure(reads, maps, cfg) {
        let backbone = group[0].0;
        for part in variant::split_on_recurrent_indels(group, backbone, reads, cfg, 0) {
            let pi = parts.len();
            parts.push(part.iter().map(|m| m.0).collect());
            for (iso, alt_start) in ends::split_by_ends(&part, reads, backbone, cfg) {
                if iso.len() >= cfg.min_iso {
                    isoforms.push((iso, alt_start, pi));
                }
            }
        }
    }
    isoforms.sort_by_key(|m| Reverse(m.0.len()));
    (isoforms, parts)
}

/// Step 1. The first member of each group is its backbone, placed at [0, length].
///
/// One round usually takes in a whole dominant structure, because every truncated copy of
/// the longest read fits it; only reads of other structures carry over. A cluster holding
/// one gene therefore costs about one test per read.
fn group_by_structure(reads: &[Read], maps: &[MinimizerMap], cfg: &Cfg) -> Vec<Vec<Member>> {
    // Longest first. Filtering keeps the order, so remaining[0] is always the next backbone.
    let mut remaining: Vec<usize> = (0..reads.len()).collect();
    remaining.sort_by_key(|&i| Reverse(reads[i].seq.len()));

    let mut groups: Vec<Vec<Member>> = Vec::new();
    while !remaining.is_empty() {
        let bb = remaining[0];
        let mut group: Vec<Member> = vec![(bb, 0, reads[bb].seq.len() as i32)];
        let mut rest: Vec<usize> = Vec::new();
        for &r in &remaining[1..] {
            match fit::fits_interval(&maps[r], &maps[bb], &reads[r].seq, &reads[bb].seq, cfg) {
                Some((start, end)) => group.push((r, start, end)),
                None => rest.push(r),
            }
        }
        groups.push(group);
        remaining = rest;
    }
    groups
}
