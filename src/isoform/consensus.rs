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
//!   * where fewer than [`MIN_OWN`] of the isoform's own reads cover a boundary or vote in a
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

use super::anchors::{chain_anchors, shared_anchors, MinimizerMap};
use super::Read;

/// Fewest reads that must carry a backbone minimizer for it to bound a window. Two reads
/// can only tie, which keeps the backbone, so a stretch covered by the backbone alone or by
/// one other read is left as it is.
const MIN_VOTERS: u32 = 2;
/// Where fewer than this many of the isoform's own reads cover a boundary or vote in a
/// window, the other reads of its structure group vote too. They share its exon structure
/// and differ only in their ends, so they carry the same sequence there. Where the isoform
/// has enough reads of its own, only they vote, so a small real difference (a shifted
/// splice site) is not outvoted by a sibling isoform.
const MIN_OWN: u32 = 3;
/// At most this many group reads help, the longest: stretches the isoform's own reads cover
/// thinly are mostly its 5' end, which only long reads reach, and placing every read of a
/// large group on every isoform's backbone would multiply the cost of this step.
const MAX_HELPERS: usize = 16;

/// Refine an isoform's consensus from its member reads (indices into `reads`),
/// reusing minimizer maps already computed for the whole cluster (`maps[ri]`).
/// `group` holds the reads of the isoform's structure group, which vote where fewer than
/// [`MIN_OWN`] members do (pass `members` to use members only).
/// Falls back to the backbone sequence for isoforms too small or too
/// minimizer-poor to vote on.
///
/// A window boundary must be a single, unambiguous position in every read that votes on it,
/// so this step only ever uses a minimizer where it occurs exactly once in the backbone and
/// once in the read; a read where the code repeats simply doesn't vote on that boundary.
pub fn refine_consensus_with_maps(
    reads: &[Read],
    members: &[usize],
    group: &[usize],
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
    let n = cand.len();
    let (depth, support) = depth_and_support(&hits, n);

    // Where members alone are too few, the rest of the structure group helps.
    let helpers: Vec<usize> = if depth.iter().any(|&d| d < MIN_OWN) {
        let mut is_member = vec![false; reads.len()];
        for &ri in members {
            is_member[ri] = true;
        }
        let mut others: Vec<usize> = group.iter().copied().filter(|&ri| !is_member[ri]).collect();
        if others.len() > MAX_HELPERS {
            others.select_nth_unstable_by_key(MAX_HELPERS - 1, |&ri| Reverse(reads[ri].seq.len()));
            others.truncate(MAX_HELPERS);
        }
        others
    } else {
        Vec::new()
    };
    let helper_hits: Vec<Vec<(u32, u32)>> =
        helpers.iter().map(|&ri| collinear_hits(&maps[ri], &maps[bb], &cand_of)).collect();
    let (h_depth, h_support) = depth_and_support(&helper_hits, n);

    let mut breaks: Vec<usize> = Vec::new(); // candidate indices, increasing
    for ci in 0..n {
        let (d, s) = if depth[ci] >= MIN_OWN {
            (depth[ci], support[ci])
        } else {
            (depth[ci] + h_depth[ci], support[ci] + h_support[ci])
        };
        let majority = d / 2 + 1; // strict majority of the local depth
        if s >= majority.max(MIN_VOTERS) {
            breaks.push(ci);
        }
    }
    if breaks.len() < 2 {
        return backbone.clone();
    }
    let m = breaks.len();
    let n_windows = m - 1;

    // Per-window vote tables (candidate segment bytes -> count): members, and the helpers
    // that join a window where fewer than MIN_OWN members vote.
    let tally = |hits: &[Vec<(u32, u32)>], ids: &[usize]| {
        let mut votes: Vec<HashMap<Vec<u8>, u32>> = vec![HashMap::new(); n_windows];
        let mut rpos: Vec<Option<u32>> = vec![None; m];
        for (h, &ri) in hits.iter().zip(ids) {
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
        votes
    };
    let mut votes = tally(&hits, members);
    if !helpers.is_empty() {
        for (own, extra) in votes.iter_mut().zip(tally(&helper_hits, &helpers)) {
            if own.values().sum::<u32>() < MIN_OWN {
                for (seg, c) in extra {
                    *own.entry(seg).or_insert(0) += c;
                }
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

/// Per candidate boundary: how many reads' chains span it, and how many hit it.
fn depth_and_support(hits: &[Vec<(u32, u32)>], n: usize) -> (Vec<u32>, Vec<u32>) {
    let mut step = vec![0i32; n + 1];
    let mut support = vec![0u32; n];
    for h in hits {
        if let (Some(&(first, _)), Some(&(last, _))) = (h.first(), h.last()) {
            step[first as usize] += 1;
            step[last as usize + 1] -= 1;
        }
        for &(ci, _) in h {
            support[ci as usize] += 1;
        }
    }
    let mut depth = Vec::with_capacity(n);
    let mut d = 0i32;
    for &x in &step[..n] {
        d += x;
        depth.push(d.max(0) as u32);
    }
    (depth, support)
}

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
        refine_consensus_with_maps(&reads, &members, &members, &maps)
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
    fn group_reads_vote_where_the_isoform_has_too_few_of_its_own() {
        // The isoform has the backbone (a 1-bp insertion at 500) and one correct full-length
        // read: a tie, so its own reads keep the error. Three reads of another isoform of the
        // same structure group cover position 500 correctly and break the tie.
        let truth = random_seq(900, 23);
        let mut backbone = truth.clone();
        backbone.insert(500, b'T');
        let seqs = vec![backbone, truth.clone(), truth[300..].to_vec(), truth[320..].to_vec(),
                        truth[350..].to_vec()];
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| Read { name: i.to_string(), seq: s.clone() })
            .collect();
        let maps: Vec<MinimizerMap> = reads.iter().map(|r| minimizer_map(&r.seq, 15, 10)).collect();
        let members = [0, 1];
        let group = [0, 1, 2, 3, 4];
        assert_ne!(refine_consensus_with_maps(&reads, &members, &members, &maps), truth);
        assert_eq!(refine_consensus_with_maps(&reads, &members, &group, &maps), truth);
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
