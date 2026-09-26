//! Step 4: consensus per isoform, from anchor-windowed modal segments.
//!
//! The reads of one isoform are near-identical and collinear (that is exactly
//! why they were grouped). So instead of a multiple-sequence alignment / POA we
//! exploit the shared unique minimizer anchors we already know how to compute:
//!
//!   * backbone      = the longest read;
//!   * each read is placed on the backbone by its collinear chain of backbone minimizers
//!     (codes that occur once in the backbone and once in the read); the chain's first and
//!     last hit give the stretch of backbone the read covers;
//!   * WINDOW BOUNDARIES are the backbone minimizers carried by a strict majority of the
//!     reads that COVER that position, and by at least [`MIN_VOTERS`] reads — an anchor a
//!     single read's error can't destroy. (Backbone-only minimizers are never boundaries:
//!     they cluster exactly at the backbone's own errors, and anchoring there would leave
//!     that error with no other read able to vote it out — the whole point is to overrule it.)
//!   * a read contributes to a window iff its chain holds both of the window's boundaries —
//!     its substring between those two anchors is that window's candidate;
//!   * the window's consensus is the MOST FREQUENT exact candidate substring, so
//!     a backbone error sits inside a majority-anchored window and is outvoted.
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

use std::collections::HashMap;

use super::anchors::{chain_anchors, shared_anchors, MinimizerMap};
use super::Read;

/// Fewest reads that must carry a backbone minimizer for it to bound a window. Two reads
/// can only tie, which keeps the backbone, so a stretch covered by the backbone alone or by
/// one other read is left as it is.
const MIN_VOTERS: u32 = 2;

/// Refine an isoform's consensus from its member reads (indices into `reads`),
/// reusing minimizer maps already computed for the whole cluster (`maps[ri]`).
/// Falls back to the backbone sequence for isoforms too small or too
/// minimizer-poor to vote on.
///
/// A window boundary must be a single, unambiguous position in every read that votes on it,
/// so this step only ever uses a minimizer where it occurs exactly once in the backbone and
/// once in the read; a read where the code repeats simply doesn't vote on that boundary.
pub fn refine_consensus_with_maps(
    reads: &[Read],
    members: &[usize],
    maps: &[MinimizerMap],
) -> Vec<u8> {
    // Backbone = longest member.
    let bb = *members.iter().max_by_key(|&&i| reads[i].seq.len()).unwrap();
    let backbone = &reads[bb].seq;
    if members.len() < 2 {
        return backbone.clone();
    }

    // Candidate boundaries: backbone minimizers unique in the backbone, by position. They
    // are k-mer starts of distinct codes, so the positions are distinct.
    let mut cand: Vec<u32> =
        maps[bb].runs().filter(|(_, p)| p.len() == 1).map(|(_, p)| p[0]).collect();
    cand.sort_unstable();
    if cand.len() < 2 {
        return backbone.clone();
    }
    let mut cand_of = vec![NONE; backbone.len()]; // backbone position -> candidate index
    for (ci, &p) in cand.iter().enumerate() {
        cand_of[p as usize] = ci as u32;
    }

    // Each member's collinear hits as (candidate index, read position), increasing in both.
    let hits: Vec<Vec<(u32, u32)>> =
        members.iter().map(|&ri| collinear_hits(&maps[ri], &maps[bb], &cand_of)).collect();

    // Local depth (members whose chain spans a candidate) and support (members hitting it).
    let n = cand.len();
    let mut depth_step = vec![0i32; n + 1];
    let mut support = vec![0u32; n];
    for h in &hits {
        if let (Some(&(first, _)), Some(&(last, _))) = (h.first(), h.last()) {
            depth_step[first as usize] += 1;
            depth_step[last as usize + 1] -= 1;
        }
        for &(ci, _) in h {
            support[ci as usize] += 1;
        }
    }
    let mut breaks: Vec<usize> = Vec::new(); // candidate indices, increasing
    let mut depth = 0i32;
    for ci in 0..n {
        depth += depth_step[ci];
        let majority = (depth.max(0) as u32) / 2 + 1; // strict majority of the local depth
        if support[ci] >= majority.max(MIN_VOTERS) {
            breaks.push(ci);
        }
    }
    if breaks.len() < 2 {
        return backbone.clone();
    }
    let m = breaks.len();
    let n_windows = m - 1;

    // Per-window vote table: candidate segment bytes -> count.
    let mut votes: Vec<HashMap<Vec<u8>, u32>> = vec![HashMap::new(); n_windows];
    let mut rpos: Vec<Option<u32>> = vec![None; m];
    for (h, &ri) in hits.iter().zip(members) {
        // Read position at each boundary, if the read's chain holds it.
        let mut k = 0;
        for (bi, &ci) in breaks.iter().enumerate() {
            while k < h.len() && (h[k].0 as usize) < ci {
                k += 1;
            }
            rpos[bi] = (k < h.len() && h[k].0 as usize == ci).then(|| h[k].1);
        }
        let seq = &reads[ri].seq;
        for wi in 0..n_windows {
            if let (Some(a), Some(b)) = (rpos[wi], rpos[wi + 1]) {
                let seg = seq[a as usize..b as usize].to_vec();
                *votes[wi].entry(seg).or_insert(0) += 1;
            }
        }
    }

    // Assemble: unanchored leading end, modal windows, unanchored trailing end.
    let at = |bi: usize| cand[breaks[bi]] as usize;
    let mut out: Vec<u8> = Vec::with_capacity(backbone.len());
    out.extend_from_slice(&backbone[..at(0)]);
    for wi in 0..n_windows {
        let default = &backbone[at(wi)..at(wi + 1)];
        out.extend_from_slice(&modal_segment(&votes[wi], default));
    }
    out.extend_from_slice(&backbone[at(m - 1)..]);
    out
}

const NONE: u32 = u32::MAX;

/// A read's collinear hits on the backbone's candidate boundaries, as (candidate index, read
/// position), increasing in both. `cand_of` maps a backbone position to its candidate index,
/// or [`NONE`]. A code that repeats in the read pairs one backbone position with several
/// read positions; such an ambiguous boundary is dropped for this read.
fn collinear_hits(read: &MinimizerMap, bb: &MinimizerMap, cand_of: &[u32]) -> Vec<(u32, u32)> {
    let mut anchors: Vec<(u32, u32)> = shared_anchors(read, bb)
        .into_iter()
        .filter(|&(_, q)| cand_of[q as usize] != NONE)
        .collect();
    anchors.sort_unstable_by_key(|&(p, q)| (q, p));
    let mut unique: Vec<(u32, u32)> = Vec::with_capacity(anchors.len());
    for run in anchors.chunk_by(|a, b| a.1 == b.1) {
        if run.len() == 1 {
            unique.push(run[0]);
        }
    }
    chain_anchors(&mut unique).into_iter().map(|(p, q)| (cand_of[q as usize], p)).collect()
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
            .map(|(i, s)| Read { name: i.to_string(), seq: s.clone() })
            .collect();
        let maps: Vec<MinimizerMap> = reads.iter().map(|r| minimizer_map(&r.seq, 15, 10)).collect();
        let members: Vec<usize> = (0..reads.len()).collect();
        refine_consensus_with_maps(&reads, &members, &maps)
    }

    #[test]
    fn corrects_a_backbone_error_that_only_a_minority_of_members_cover() {
        // A 900-bp transcript. The backbone (the longest read) has a 1-bp insertion at 60;
        // two more full-length reads are correct there; six reads are 5'-truncated and start
        // at 400, so only 3 of 9 members cover position 60. A majority of all members (5) is
        // out of reach there, but 2 of the 3 reads that cover it agree.
        let truth = random_seq(900, 7);
        let mut backbone = truth.clone();
        backbone.insert(60, b'A');
        let mut seqs = vec![backbone, truth.clone(), truth[5..].to_vec()];
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
    fn a_truncated_read_error_does_not_enter_the_consensus() {
        // An error in one truncated read is outvoted by the reads covering the same window.
        let truth = random_seq(700, 3);
        let mut bad = truth[200..].to_vec();
        bad.insert(150, b'C');
        let seqs = vec![truth.clone(), truth[100..].to_vec(), truth[150..].to_vec(), bad];
        assert_eq!(consensus_of(&seqs), truth);
    }
}
