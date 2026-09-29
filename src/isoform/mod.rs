//! Isoform detection: split the reads of one cluster into isoforms, each with a consensus
//! sequence. Used by `predict` after clustering, and by `find-isoforms` on an existing
//! cluster assignment.
//!
//! [`detect_isoforms`] runs four steps on a cluster:
//!
//!   1. Group by structure ([`fit`]). Reads are taken longest first: each round the longest
//!      unassigned read is the backbone, every other unassigned read that fits it joins its
//!      group, and the rest wait for a later round.
//!   2. Split each group by where its reads start and end ([`ends`]).
//!   3. Keep groups of at least `min_iso` reads and build each one's consensus ([`consensus`]).
//!   4. Merge isoforms whose consensuses show they are the same transcript ([`collapse`]).
//!
//! Which test tells two isoforms apart:
//!
//!   extra sequence in one piece (exon, retained intron)   interior-gap indel run   step 1
//!   same-length but different sequence (MXE)              interior-gap identity    step 1
//!   different sequence at one end                         unmatched flank          step 1
//!   different transcription start or end                  end peaks                step 2
//!     (with `--start-split off|auto` 3' ends only, plus validated starts for auto: 5'
//!     starts in cDNA mostly mark truncation)
//!   the same transcript found by two groups               consensus identity       step 4
//!
//! Differences confined to a flank shorter than `max_flank` are not detected. The defaults
//! were tuned on PacBio HiFi SIRV reads; ONT data has not been measured.

use std::cmp::Reverse;
use std::time::Instant;

mod align; // block-aligner calls: flank extension, global alignment, longest indel
mod anchors; // minimizer anchors and their collinear chain
mod collapse; // step 4
mod consensus; // step 3
mod ends; // step 2
mod find_isoforms; // the `find-isoforms` subcommand
mod fit; // step 1
mod options; // IsoformOptions, its defaults, flags and help text
mod output; // output files
pub(crate) mod stats; // diagnostic counters

use anchors::MinimizerMap;
pub use find_isoforms::{parse_args, run, print_usage, Config};
pub use options::{IsoformOptions, StartSplit};
pub(crate) use options::{apply_isoform_flag, ISOFORM_HELP};
pub(crate) use output::{format_cluster_records, OutputFiles, ClusterRecords};

pub struct Read {
    pub name: String,
    pub seq: Vec<u8>,
}

/// One detected isoform.
pub(crate) struct Isoform {
    /// The isoform's reads, as indices into the cluster's reads.
    pub members: Vec<usize>,
    /// Consensus sequence.
    pub consensus: Vec<u8>,
    /// Starts at a downstream transcription start that `--start-split auto` validated, so
    /// collapse must not fold it into a longer isoform that differs only at the 5' end.
    pub alt_start: bool,
}

/// A read placed on its group's backbone: (read_idx, start, end), with start and end in
/// backbone coordinates.
type PlacedRead = (usize, i32, i32);

/// Detect the isoforms in one cluster's reads. Clusters are independent, so callers can
/// detect many at once.
pub(crate) fn detect_isoforms(reads: &[Read], options: &IsoformOptions) -> Vec<Isoform> {
    let lap_start = Instant::now();
    let read_maps: Vec<MinimizerMap> = reads
        .iter()
        .map(|read| anchors::minimizer_map(&read.seq, options.minimizer_k, options.minimizer_w))
        .collect();
    stats::add_elapsed(&stats::T_MAPS, lap_start);

    let lap_start = Instant::now();
    let (drafts, reads_by_group) = partition(reads, &read_maps, options);
    stats::add_elapsed(&stats::T_PARTITION, lap_start);

    let lap_start = Instant::now();
    let isoforms: Vec<Isoform> = drafts
        .into_iter()
        .map(|(members, alt_start, group_index)| {
            let consensus = consensus::build_consensus(
                reads,
                &members,
                &reads_by_group[group_index],
                &read_maps,
            );
            Isoform { members, consensus, alt_start }
        })
        .collect();
    stats::add_elapsed(&stats::T_CONS, lap_start);

    collapse::merge_duplicates(isoforms, reads, &read_maps, options)
}

/// Steps 1–2 and the `min_iso` filter. Returns the drafts (isoforms that have reads but no
/// consensus yet), largest first, each as its read indices, whether it starts at a validated
/// downstream start, and the `group_index` of its structure group in the second list, which
/// holds each group's reads: they share the isoform's structure, so they may vote in its
/// consensus.
fn partition(
    reads: &[Read],
    read_maps: &[MinimizerMap],
    options: &IsoformOptions,
) -> (Vec<(Vec<usize>, bool, usize)>, Vec<Vec<usize>>) {
    let mut drafts: Vec<(Vec<usize>, bool, usize)> = Vec::new();
    let mut reads_by_group: Vec<Vec<usize>> = Vec::new();
    for group in &group_by_structure(reads, read_maps, options) {
        let (backbone_idx, _, _) = group[0];
        let group_index = reads_by_group.len();
        reads_by_group.push(group.iter().map(|&(read_idx, _, _)| read_idx).collect());
        for (isoform, alt_start) in ends::split_by_ends(group, reads, backbone_idx, options) {
            if isoform.len() >= options.min_iso {
                drafts.push((isoform, alt_start, group_index));
            }
        }
    }
    drafts.sort_by_key(|(members, _, _)| Reverse(members.len()));
    (drafts, reads_by_group)
}

/// Step 1. The first read of each group is its backbone, placed at [0, length].
///
/// One round usually takes in a whole dominant structure, because every truncated copy of
/// the longest read fits it; only reads of other structures carry over. A cluster holding
/// one gene therefore costs about one test per read.
fn group_by_structure(
    reads: &[Read],
    read_maps: &[MinimizerMap],
    options: &IsoformOptions,
) -> Vec<Vec<PlacedRead>> {
    // Longest first. Filtering keeps the order, so remaining[0] is always the next backbone.
    let mut remaining: Vec<usize> = (0..reads.len()).collect();
    remaining.sort_by_key(|&read_idx| Reverse(reads[read_idx].seq.len()));

    let mut groups: Vec<Vec<PlacedRead>> = Vec::new();
    while !remaining.is_empty() {
        let backbone_idx = remaining[0];
        let mut group: Vec<PlacedRead> =
            vec![(backbone_idx, 0, reads[backbone_idx].seq.len() as i32)];
        let mut rest: Vec<usize> = Vec::new();
        for &read_idx in &remaining[1..] {
            match fit::place_on_backbone(
                &read_maps[read_idx],
                &read_maps[backbone_idx],
                &reads[read_idx].seq,
                &reads[backbone_idx].seq,
                options,
            ) {
                Some((start, end)) => group.push((read_idx, start, end)),
                None => rest.push(read_idx),
            }
        }
        groups.push(group);
        remaining = rest;
    }
    groups
}
