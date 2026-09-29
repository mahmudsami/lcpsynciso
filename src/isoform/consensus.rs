//! Step 3: consensus per isoform, from anchor-windowed modal segments.
//!
//! The reads of one isoform are near-identical and collinear (that is exactly
//! why they were grouped). So instead of a multiple-sequence alignment / POA we
//! exploit the shared unique minimizer anchors we already know how to compute:
//!
//!   * backbone      = the longest read;
//!   * each read is placed on the backbone by its collinear chain of backbone minimizers
//!     (k-mers that occur once in the backbone and once in the read); the chain's first and
//!     last hit give the stretch of backbone the read covers;
//!   * WINDOW BOUNDARIES are the backbone minimizers carried by a strict majority of the
//!     reads that COVER that position, and by at least [`MIN_BOUNDARY_SUPPORT`] reads — an
//!     anchor a single read's error can't destroy. (Backbone-only minimizers are never
//!     boundaries: they cluster exactly at the backbone's own errors, and anchoring there
//!     would leave that error with no other read able to vote it out — the whole point is to
//!     overrule it.)
//!   * a read contributes to a window iff its chain holds both of the window's boundaries —
//!     its substring between those two anchors is that window's candidate;
//!   * the window's consensus is the MOST FREQUENT exact candidate substring, so
//!     a backbone error sits inside a majority-anchored window and is outvoted.
//!   * where fewer than [`MIN_OWN_READS`] of the isoform's own reads cover a boundary or vote in a
//!     window, the other reads of its structure group vote too: same exons, other ends, so
//!     the same sequence there (e.g. the full-length isoform borrowing the reads of an
//!     isoform at a downstream start).
//!
//! The majority is local because most reads of an isoform are usually 5'-truncated. With a
//! majority over ALL members, the 5' stretch that only the few full-length reads cover had
//! no boundaries at all and was copied verbatim from the backbone, errors included. On UHRR
//! HiFi that was the main reason an isoform's coding sequence carried a 1-bp frameshift.
//!
//! Because the flanking k-mers are exact matches in every contributing read, the
//! segments are perfectly registered and concatenate seamlessly — no inner
//! alignment at all. Cost is ~O(total read length): one minimizer merge per read
//! plus hashing.
//!
//! Why the mode is the truth: at window length L and per-base error e, a read is
//! error-free across a window with probability (1-e)^L — for HiFi (e~0.2%) and
//! L~w=10 that is ~0.98, so the overwhelming majority of reads vote the true
//! string. This also fixes the dominant HiFi error mode (homopolymer indels),
//! since an indel changes the segment string and stays a minority. Reads that
//! are truncated or lost an anchor to a local error simply skip that window; the
//! rest still carry it.

use std::cmp::Reverse;
use std::collections::HashMap;

use super::anchors::{chain_anchors, shared_anchors, Anchor, MinimizerMap};
use super::Read;

/// Fewest reads that must carry a backbone minimizer for it to bound a window. Two reads
/// can only tie, which keeps the backbone, so a stretch covered by the backbone alone or by
/// one other read is left as it is.
const MIN_BOUNDARY_SUPPORT: u32 = 2;
/// Where fewer than this many of the isoform's own reads cover a boundary or vote in a
/// window, the other reads of its structure group vote too. They share its exon structure
/// and differ only in their ends, so they carry the same sequence there. Where the isoform
/// has enough reads of its own, only they vote, so a small real difference (a shifted
/// splice site) is not outvoted by a sibling isoform.
const MIN_OWN_READS: u32 = 3;
/// At most this many group reads help, the longest: stretches the isoform's own reads cover
/// thinly are mostly its 5' end, which only long reads reach, and placing every read of a
/// large group on every isoform's backbone would multiply the cost of this step.
const MAX_HELPERS: usize = 16;

/// Build an isoform's consensus from its member reads (indices into `reads`),
/// reusing minimizer maps already computed for the whole cluster (`read_maps[read_idx]`).
/// `group_reads` holds the reads of the isoform's structure group, which vote where fewer than
/// [`MIN_OWN_READS`] members do (pass `members` to use members only).
/// Falls back to the backbone sequence for isoforms too small or too
/// minimizer-poor to vote on.
///
/// A window boundary must be a single, unambiguous position in every read that votes on it,
/// so this step only ever uses a minimizer where it occurs exactly once in the backbone and
/// once in the read; a read where the k-mer repeats simply doesn't vote on that boundary.
pub fn build_consensus(
    reads: &[Read],
    members: &[usize],
    group_reads: &[usize],
    read_maps: &[MinimizerMap],
) -> Vec<u8> {
    // Backbone = longest member.
    let backbone_idx = *members.iter().max_by_key(|&&read_idx| reads[read_idx].seq.len()).unwrap();
    let backbone_seq = &reads[backbone_idx].seq;
    if members.len() < 2 {
        return backbone_seq.clone();
    }

    // Candidate boundaries: backbone minimizers unique in the backbone, by position. They
    // are k-mer starts of distinct k-mers, so the positions are distinct.
    let mut candidates: Vec<u32> =
        read_maps[backbone_idx].runs().filter(|(_, p)| p.len() == 1).map(|(_, p)| p[0]).collect();
    candidates.sort_unstable();
    if candidates.len() < 2 {
        return backbone_seq.clone();
    }
    // Backbone position -> candidate index.
    let mut candidate_at = vec![NO_CANDIDATE; backbone_seq.len()];
    for (ci, &backbone_pos) in candidates.iter().enumerate() {
        candidate_at[backbone_pos as usize] = ci as u32;
    }

    // Each member's collinear hits as (candidate index, read position), increasing in both.
    let hits: Vec<Vec<(u32, u32)>> = members
        .iter()
        .map(|&read_idx| {
            collinear_hits(&read_maps[read_idx], &read_maps[backbone_idx], &candidate_at)
        })
        .collect();
    let n_candidates = candidates.len();
    let (depth, support) = depth_and_support(&hits, n_candidates);

    // Where members alone are too few, the rest of the structure group helps.
    let helpers: Vec<usize> = if depth.iter().any(|&d| d < MIN_OWN_READS) {
        let mut is_member = vec![false; reads.len()];
        for &read_idx in members {
            is_member[read_idx] = true;
        }
        let mut others: Vec<usize> =
            group_reads.iter().copied().filter(|&read_idx| !is_member[read_idx]).collect();
        if others.len() > MAX_HELPERS {
            others.select_nth_unstable_by_key(MAX_HELPERS - 1, |&read_idx| {
                Reverse(reads[read_idx].seq.len())
            });
            others.truncate(MAX_HELPERS);
        }
        others
    } else {
        Vec::new()
    };
    let helper_hits: Vec<Vec<(u32, u32)>> = helpers
        .iter()
        .map(|&read_idx| {
            collinear_hits(&read_maps[read_idx], &read_maps[backbone_idx], &candidate_at)
        })
        .collect();
    let (helper_depth, helper_support) = depth_and_support(&helper_hits, n_candidates);

    let mut boundaries: Vec<usize> = Vec::new(); // candidate indices, increasing
    for ci in 0..n_candidates {
        let (local_depth, local_support) = if depth[ci] >= MIN_OWN_READS {
            (depth[ci], support[ci])
        } else {
            (depth[ci] + helper_depth[ci], support[ci] + helper_support[ci])
        };
        let majority = local_depth / 2 + 1; // strict majority of the local depth
        if local_support >= majority.max(MIN_BOUNDARY_SUPPORT) {
            boundaries.push(ci);
        }
    }
    if boundaries.len() < 2 {
        return backbone_seq.clone();
    }
    let n_boundaries = boundaries.len();
    let n_windows = n_boundaries - 1;

    // Per-window vote tables (candidate segment bytes -> count): members, and the helpers
    // that join a window where fewer than MIN_OWN_READS members vote.
    let tally = |hits: &[Vec<(u32, u32)>], read_indices: &[usize]| {
        let mut votes: Vec<HashMap<Vec<u8>, u32>> = vec![HashMap::new(); n_windows];
        let mut read_pos_at: Vec<Option<u32>> = vec![None; n_boundaries];
        for (h, &read_idx) in hits.iter().zip(read_indices) {
            // Read position at each boundary, if the read's chain holds it.
            let mut hit_idx = 0;
            for (bi, &ci) in boundaries.iter().enumerate() {
                while hit_idx < h.len() && (h[hit_idx].0 as usize) < ci {
                    hit_idx += 1;
                }
                read_pos_at[bi] =
                    (hit_idx < h.len() && h[hit_idx].0 as usize == ci).then(|| h[hit_idx].1);
            }
            let seq = &reads[read_idx].seq;
            for wi in 0..n_windows {
                if let (Some(read_start), Some(read_end)) = (read_pos_at[wi], read_pos_at[wi + 1]) {
                    let seg = seq[read_start as usize..read_end as usize].to_vec();
                    *votes[wi].entry(seg).or_insert(0) += 1;
                }
            }
        }
        votes
    };
    let mut votes = tally(&hits, members);
    if !helpers.is_empty() {
        for (own, extra) in votes.iter_mut().zip(tally(&helper_hits, &helpers)) {
            if own.values().sum::<u32>() < MIN_OWN_READS {
                for (seg, c) in extra {
                    *own.entry(seg).or_insert(0) += c;
                }
            }
        }
    }

    // Assemble: unanchored leading end, modal windows, unanchored trailing end.
    let boundary_pos = |bi: usize| candidates[boundaries[bi]] as usize;
    let mut out: Vec<u8> = Vec::with_capacity(backbone_seq.len());
    out.extend_from_slice(&backbone_seq[..boundary_pos(0)]);
    for wi in 0..n_windows {
        let default = &backbone_seq[boundary_pos(wi)..boundary_pos(wi + 1)];
        out.extend_from_slice(&modal_segment(&votes[wi], default));
    }
    out.extend_from_slice(&backbone_seq[boundary_pos(n_boundaries - 1)..]);
    out
}

const NO_CANDIDATE: u32 = u32::MAX;

/// Per candidate boundary: how many reads' chains span it, and how many hit it.
fn depth_and_support(hits: &[Vec<(u32, u32)>], n_candidates: usize) -> (Vec<u32>, Vec<u32>) {
    let mut step = vec![0i32; n_candidates + 1];
    let mut support = vec![0u32; n_candidates];
    for h in hits {
        if let (Some(&(first, _)), Some(&(last, _))) = (h.first(), h.last()) {
            step[first as usize] += 1;
            step[last as usize + 1] -= 1;
        }
        for &(ci, _) in h {
            support[ci as usize] += 1;
        }
    }
    let mut depth = Vec::with_capacity(n_candidates);
    let mut d = 0i32;
    for &x in &step[..n_candidates] {
        d += x;
        depth.push(d.max(0) as u32);
    }
    (depth, support)
}

/// A read's collinear hits on the backbone's candidate boundaries, as (candidate index, read
/// position), increasing in both. `candidate_at` maps a backbone position to its candidate
/// index, or [`NO_CANDIDATE`]. A k-mer that repeats in the read pairs one backbone position with
/// several read positions; such an ambiguous boundary is dropped for this read.
fn collinear_hits(
    read_map: &MinimizerMap,
    backbone_map: &MinimizerMap,
    candidate_at: &[u32],
) -> Vec<(u32, u32)> {
    let mut anchors: Vec<Anchor> = shared_anchors(read_map, backbone_map)
        .into_iter()
        .filter(|&(_, backbone_pos)| candidate_at[backbone_pos as usize] != NO_CANDIDATE)
        .collect();
    anchors.sort_unstable_by_key(|&(read_pos, backbone_pos)| (backbone_pos, read_pos));
    let mut unique: Vec<Anchor> = Vec::with_capacity(anchors.len());
    for run in anchors.chunk_by(|a, b| a.1 == b.1) {
        if run.len() == 1 {
            unique.push(run[0]);
        }
    }
    chain_anchors(&mut unique)
        .into_iter()
        .map(|(read_pos, backbone_pos)| (candidate_at[backbone_pos as usize], read_pos))
        .collect()
}

/// Most frequent segment in a window; ties break toward the backbone default
/// (conservative — never change the backbone without a strict majority), then
/// toward the lexicographically smaller bytes for determinism. An empty table
/// (no read spanned the window) yields the backbone default.
fn modal_segment(table: &HashMap<Vec<u8>, u32>, default: &[u8]) -> Vec<u8> {
    let mut best: Option<(&[u8], u32)> = None;
    for (seg, &c) in table {
        let cand = seg.as_slice();
        best = Some(match best {
            None => (cand, c),
            Some((bseg, bc)) => {
                if c > bc {
                    (cand, c)
                } else if c < bc {
                    (bseg, bc)
                } else if bseg == default {
                    (bseg, bc)
                } else if cand == default || cand < bseg {
                    (cand, c)
                } else {
                    (bseg, bc)
                }
            }
        });
    }
    match best {
        Some((seg, _)) => seg.to_vec(),
        None => default.to_vec(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::isoform::anchors::minimizer_map;

    /// A fixed pseudo-random sequence, so minimizers are unique and well spread.
    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    fn consensus_of(seqs: &[Vec<u8>]) -> Vec<u8> {
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(read_idx, seq)| Read { name: read_idx.to_string(), seq: seq.clone() })
            .collect();
        let read_maps: Vec<MinimizerMap> =
            reads.iter().map(|r| minimizer_map(&r.seq, 15, 10)).collect();
        let members: Vec<usize> = (0..reads.len()).collect();
        build_consensus(&reads, &members, &members, &read_maps)
    }

    #[test]
    fn corrects_a_backbone_error_that_only_a_minority_of_members_cover() {
        // A 900-bp transcript. The backbone (the longest read) has a 1-bp insertion at 60;
        // two more full-length reads are correct there; six reads are 5'-truncated and start
        // at 400, so only 3 of 9 members cover position 60. A majority of all members (5) is
        // out of reach there, but 2 of the 3 reads that cover it agree.
        let truth = random_seq(900, 7);
        let mut backbone_seq = truth.clone();
        backbone_seq.insert(60, b'A');
        let mut seqs = vec![backbone_seq, truth.clone(), truth[5..].to_vec()];
        for _ in 0..6 {
            seqs.push(truth[400..].to_vec());
        }
        assert_eq!(consensus_of(&seqs), truth);
    }

    #[test]
    fn keeps_the_backbone_where_only_one_other_read_disagrees() {
        // Two reads that disagree have no majority: the backbone is kept.
        let truth = random_seq(600, 11);
        let mut other = truth.clone();
        other.remove(300);
        assert_eq!(consensus_of(&[truth.clone(), other]), truth);
    }

    #[test]
    fn group_reads_vote_where_the_isoform_has_too_few_of_its_own() {
        // The isoform has the backbone (a 1-bp insertion at 500) and one correct full-length
        // read: a tie, so its own reads keep the error. Three reads of another isoform of the
        // same structure group cover position 500 correctly and break the tie.
        let truth = random_seq(900, 23);
        let mut backbone_seq = truth.clone();
        backbone_seq.insert(500, b'T');
        let seqs = vec![backbone_seq, truth.clone(), truth[300..].to_vec(), truth[320..].to_vec(),
                        truth[350..].to_vec()];
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(read_idx, seq)| Read { name: read_idx.to_string(), seq: seq.clone() })
            .collect();
        let read_maps: Vec<MinimizerMap> =
            reads.iter().map(|r| minimizer_map(&r.seq, 15, 10)).collect();
        let members = [0, 1];
        let group_reads = [0, 1, 2, 3, 4];
        assert_ne!(build_consensus(&reads, &members, &members, &read_maps), truth);
        assert_eq!(build_consensus(&reads, &members, &group_reads, &read_maps), truth);
    }

    #[test]
    fn a_truncated_read_error_does_not_enter_the_consensus() {
        // An error in one truncated read is outvoted by the reads covering the same window.
        let truth = random_seq(700, 3);
        let mut bad = truth[200..].to_vec();
        bad.insert(150, b'C');
        let seqs = vec![truth.clone(), truth[100..].to_vec(), truth[150..].to_vec(), bad];
        assert_eq!(consensus_of(&seqs), truth);
    }
}
