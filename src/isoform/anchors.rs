//! Minimizer anchors between two sequences, and the collinear chain through them.
//!
//! An anchor pairs a read position with a backbone position that share a minimizer k-mer.
//! Every occurrence of a minimizer is kept, not only k-mers unique to one sequence, so a
//! minimizer repeated in a read or backbone can seed more than one candidate anchor; the
//! chain step below is what discards the pairings that aren't part of a real collinear
//! alignment. The chain is the largest set of anchors whose positions increase together in
//! both sequences.

use std::cell::RefCell;
use std::cmp::{Ordering, Reverse};

/// A shared minimizer pair, as `(read_pos, backbone_pos)`.
pub(super) type Anchor = (u32, u32);

// ── Minimizers ───────────────────────────────────────────────────────────────

/// A sequence's minimizers: every minimizer k-mer with every position where it is the window
/// minimum, as two parallel arrays sorted by k-mer, then position. A k-mer is packed into a
/// `u64`, two bits per base. About 12 bytes per minimizer, in
/// two allocations; detection holds one of these for every read of every cluster in flight.
#[derive(Debug)]
pub(super) struct MinimizerMap {
    kmers: Box<[u64]>,
    positions: Box<[u32]>,
}

impl MinimizerMap {
    /// Positions of `kmer`, increasing; empty if the sequence lacks it.
    #[cfg(test)]
    pub(super) fn get(&self, kmer: u64) -> &[u32] {
        let lo = self.kmers.partition_point(|&k| k < kmer);
        let hi = lo + self.kmers[lo..].partition_point(|&k| k == kmer);
        &self.positions[lo..hi]
    }

    /// Each distinct k-mer with its positions, in k-mer order.
    pub(super) fn runs(&self) -> impl Iterator<Item = (u64, &[u32])> {
        let mut i = 0;
        std::iter::from_fn(move || {
            let kmer = *self.kmers.get(i)?;
            let end = run_end(&self.kmers, i);
            let run = (kmer, &self.positions[i..end]);
            i = end;
            Some(run)
        })
    }
}

/// End of the run of equal k-mers starting at `i`.
fn run_end(kmers: &[u64], i: usize) -> usize {
    let mut end = i + 1;
    while end < kmers.len() && kmers[end] == kmers[i] {
        end += 1;
    }
    end
}

/// Every (read position, backbone position) pair at which the two sequences share a
/// minimizer k-mer: each shared k-mer pairs each of its read positions with each of its
/// backbone positions. A merge of the two sorted k-mer arrays, so no hashing.
pub(super) fn shared_anchors(read_map: &MinimizerMap, backbone_map: &MinimizerMap) -> Vec<Anchor> {
    let (read_kmers, backbone_kmers) = (&read_map.kmers, &backbone_map.kmers);
    let (mut i, mut j) = (0, 0);
    let mut anchors = Vec::new();
    while i < read_kmers.len() && j < backbone_kmers.len() {
        match read_kmers[i].cmp(&backbone_kmers[j]) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let (i_end, j_end) = (run_end(read_kmers, i), run_end(backbone_kmers, j));
                for &read_pos in &read_map.positions[i..i_end] {
                    for &backbone_pos in &backbone_map.positions[j..j_end] {
                        anchors.push((read_pos, backbone_pos));
                    }
                }
                i = i_end;
                j = j_end;
            }
        }
    }
    anchors
}

/// Packed k-mer (two bits per base) and start position of every k-mer; k-mers containing a
/// non-ACGT base are skipped.
fn encode_kmers(seq: &[u8], k: usize) -> Vec<(u64, u32)> {
    let mut out = Vec::with_capacity(seq.len());
    let mask: u64 = if k >= 32 { u64::MAX } else { (1u64 << (2 * k)) - 1 };
    let (mut val, mut len) = (0u64, 0usize);
    for (i, &b) in seq.iter().enumerate() {
        let c = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => {
                len = 0;
                val = 0;
                continue;
            }
        };
        val = ((val << 2) | c) & mask;
        len += 1;
        if len >= k {
            out.push((val, (i + 1 - k) as u32));
        }
    }
    out
}

/// Window-`w` minimizers of `seq`: each k-mer with every position where it is the window
/// minimum. Most k-mers carry exactly one position; a k-mer the window selects more than once
/// (a repeated minimizer) carries all of them.
pub(super) fn minimizer_map(seq: &[u8], k: usize, w: usize) -> MinimizerMap {
    let all_kmers = encode_kmers(seq, k);
    let mut minimizers: Vec<(u64, u32)> = Vec::new();
    if all_kmers.len() < w.max(1) {
        minimizers = all_kmers;
    } else {
        let (mut last_kmer, mut last_pos) = (u64::MAX, u32::MAX);
        for window in all_kmers.windows(w) {
            let window_min =
                *window.iter().min_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1))).unwrap();
            if window_min.0 != last_kmer || window_min.1 != last_pos {
                minimizers.push(window_min);
                last_kmer = window_min.0;
                last_pos = window_min.1;
            }
        }
    }
    minimizers.retain(|&(kmer, _)| !is_homopolymer(kmer, k));
    minimizers.sort_unstable();
    MinimizerMap {
        kmers: minimizers.iter().map(|&(kmer, _)| kmer).collect(),
        positions: minimizers.iter().map(|&(_, pos)| pos).collect(),
    }
}

/// Whether every base of `kmer` is the same, as in a polyA tail. Such a k-mer is the
/// minimizer at every position of its run, so keeping all of them would pair each position of
/// one run with each of the other's: a ladder of equally valid anchors that the chain follows
/// to an arbitrary offset, putting the terminal anchor inside the tail instead of at the last
/// real sequence.
fn is_homopolymer(kmer: u64, k: usize) -> bool {
    let base = kmer & 3;
    (0..k).all(|i| (kmer >> (2 * i)) & 3 == base)
}

// ── Chaining ─────────────────────────────────────────────────────────────────

/// The longest chain through `anchors` whose read and backbone positions both strictly
/// increase together. Sorts `anchors` in place: by read position ascending, then by backbone
/// position descending among ties, so that anchors sharing a read position — a minimizer
/// repeated in the read — can contribute at most one to the chain (the standard ordering for
/// a longest chain strictly increasing in two coordinates, as opposed to a plain LIS).
///
/// Maximising the number of anchors keeps a mispaired anchor out of the chain, because
/// keeping it would block every anchor it sits in front of. Most comparisons have no
/// mispaired anchor, so the greedy scan runs first: when it keeps every anchor, the full set
/// is already the longest chain.
pub(super) fn chain_anchors(anchors: &mut [Anchor]) -> Vec<Anchor> {
    anchors.sort_unstable_by_key(|&(read_pos, backbone_pos)| (read_pos, Reverse(backbone_pos)));
    let chain = greedy_chain(anchors);
    if chain.len() == anchors.len() {
        return chain;
    }
    lis_chain(anchors)
}

/// Keep each anchor whose backbone position is above the last one kept. Anchors that share a
/// read position are sorted backbone-descending, so at most the first one seen can pass.
fn greedy_chain(anchors: &[Anchor]) -> Vec<Anchor> {
    let mut chain: Vec<Anchor> = Vec::with_capacity(anchors.len());
    for &(read_pos, backbone_pos) in anchors {
        match chain.last() {
            Some(&(_, last_backbone_pos)) if backbone_pos <= last_backbone_pos => continue,
            _ => chain.push((read_pos, backbone_pos)),
        }
    }
    chain
}

/// Index buffers for [`lis_chain`], reused across calls on each worker thread because
/// `place_on_backbone` chains millions of anchor sets.
#[derive(Default)]
struct ChainScratch {
    /// `tails[i]`: among chains of length `i + 1` so far, the anchor that ends one at the
    /// lowest backbone position.
    tails: Vec<u32>,
    /// Predecessor of each anchor in its best chain, or [`NO_PREV`].
    prev: Vec<u32>,
}

thread_local! {
    static SCRATCH: RefCell<ChainScratch> = RefCell::new(ChainScratch::default());
}

const NO_PREV: u32 = u32::MAX;

/// Longest strictly increasing subsequence by backbone position, in O(n log n) (patience
/// sorting with predecessor links). Relies on the ordering [`chain_anchors`] sorts into
/// `anchors` (read position ascending, backbone position descending among ties): a chain
/// strictly increasing in backbone position is then automatically strictly increasing in
/// read position too, so anchors sharing a read or backbone position never both appear in it.
fn lis_chain(anchors: &[Anchor]) -> Vec<Anchor> {
    if anchors.is_empty() {
        return Vec::new();
    }
    SCRATCH.with(|scratch| {
        let scratch = &mut *scratch.borrow_mut();
        scratch.tails.clear();
        scratch.prev.clear();
        scratch.prev.resize(anchors.len(), NO_PREV);

        for (idx, &(_, backbone_pos)) in anchors.iter().enumerate() {
            // Strictly increasing: the first tail whose backbone position is >= backbone_pos.
            let pos = scratch.tails.partition_point(|&t| anchors[t as usize].1 < backbone_pos);
            if pos > 0 {
                scratch.prev[idx] = scratch.tails[pos - 1];
            }
            if pos == scratch.tails.len() {
                scratch.tails.push(idx as u32);
            } else {
                scratch.tails[pos] = idx as u32;
            }
        }
        let mut out = Vec::with_capacity(scratch.tails.len());
        let mut cur = *scratch.tails.last().unwrap();
        loop {
            out.push(anchors[cur as usize]);
            if scratch.prev[cur as usize] == NO_PREV {
                break;
            }
            cur = scratch.prev[cur as usize];
        }
        out.reverse();
        out
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn minimizer_map_keeps_every_occurrence_of_a_repeated_kmer() {
        // "AAAC" (packed value 1) is the window minimum at both position 0 and position 10. The
        // old unique-only filter dropped such a k-mer entirely; every occurrence is kept now.
        let map = minimizer_map(b"AAACTTTTTTAAACTTTTTT", 4, 2);
        assert_eq!(map.get(1), &[0u32, 10][..]);
    }

    #[test]
    fn shared_anchors_pairs_every_occurrence_of_each_shared_kmer() {
        // K-mer 1 at read 0 and 10 and backbone 5; k-mer 7 at read 3 and backbone 8 and 20;
        // k-mer 9 only in the read and k-mer 4 only in the backbone.
        let read_map = MinimizerMap {
            kmers: vec![1, 1, 7, 9].into(),
            positions: vec![0, 10, 3, 6].into(),
        };
        let backbone_map = MinimizerMap {
            kmers: vec![1, 4, 7, 7].into(),
            positions: vec![5, 2, 8, 20].into(),
        };
        let mut anchors = shared_anchors(&read_map, &backbone_map);
        anchors.sort_unstable();
        assert_eq!(anchors, vec![(0, 5), (3, 8), (3, 20), (10, 5)]);
        assert_eq!(
            read_map.runs().map(|(kmer, positions)| (kmer, positions.len())).collect::<Vec<_>>(),
            vec![(1, 2), (7, 1), (9, 1)]
        );
    }

    #[test]
    fn minimizer_map_excludes_homopolymers() {
        // Every k-mer of a polyA run is the same, one per position: a ladder of equally
        // valid anchors that drags a chain to an arbitrary offset inside the tail.
        let map = minimizer_map(b"AAAAAAAA", 3, 4);
        assert!(map.runs().next().is_none(), "the all-A 3-mer must not be an anchor, got {map:?}");
        // The polyT 4-mer in the sequence above is excluded for the same reason.
        let map = minimizer_map(b"AAACTTTTTTAAACTTTTTT", 4, 2);
        assert!(map.get(255).is_empty(), "polyT must not be an anchor");
    }

    #[test]
    fn chain_anchors_plain_lis_unchanged() {
        // No repeated read or backbone position: behaves exactly like the old unique-only
        // chain. (2, 5) is a mispaired anchor: keeping it would block every anchor after it,
        // so the unique longest chain drops it instead.
        let mut anchors = vec![(0, 10), (1, 20), (2, 5), (3, 30), (4, 40)];
        let chain = chain_anchors(&mut anchors);
        assert_eq!(chain, vec![(0, 10), (1, 20), (3, 30), (4, 40)]);
    }

    #[test]
    fn chain_anchors_drops_extra_anchors_from_a_repeated_read_minimizer() {
        // Read position 10 pairs with two backbone candidates (10 -> 100, 10 -> 200): a
        // minimizer repeated in the backbone. At most one can survive collinear chaining.
        let mut anchors = vec![(10, 100), (10, 200), (20, 300)];
        let chain = chain_anchors(&mut anchors);
        assert_eq!(chain.len(), 2);
        // Every read position in the chain is distinct and both coordinates strictly
        // increase.
        assert!(chain.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1));
        assert_eq!(chain.last(), Some(&(20, 300)));
    }

    #[test]
    fn chain_anchors_handles_ties_on_both_sides() {
        // A minimizer repeated in the read (read_pos=30 vs two backbone hits) AND one repeated
        // in the backbone (backbone_pos=100 vs two read hits) in the same anchor set.
        let mut anchors = vec![(10, 100), (30, 110), (30, 90), (50, 100)];
        let chain = chain_anchors(&mut anchors);
        assert_eq!(chain.len(), 2, "no length-3 collinear chain exists among these anchors");
        assert!(
            chain.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1),
            "chain must strictly increase in both coordinates: {chain:?}"
        );
    }
}
