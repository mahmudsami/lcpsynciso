//! Minimizer anchors between two sequences, and the collinear chain through them.
//!
//! An anchor pairs a read position with a backbone position that share a minimizer code.
//! Every occurrence of a minimizer is kept, not only codes unique to one sequence, so a
//! minimizer repeated in a read or backbone can seed more than one candidate anchor; the
//! chain step below is what discards the pairings that aren't part of a real collinear
//! alignment. The chain is the largest set of anchors whose positions increase together in
//! both sequences.

use std::cell::RefCell;
use std::cmp::{Ordering, Reverse};

// ── Minimizers ───────────────────────────────────────────────────────────────

/// A sequence's minimizers: every code with every position where it is the window minimum,
/// as two parallel arrays sorted by code, then position. About 12 bytes per minimizer, in
/// two allocations; resolution holds one of these for every read of every cluster in flight.
#[derive(Debug)]
pub(super) struct MinimizerMap {
    codes: Box<[u64]>,
    pos: Box<[u32]>,
}

impl MinimizerMap {
    /// Positions of `code`, increasing; empty if the sequence lacks it.
    pub(super) fn get(&self, code: u64) -> &[u32] {
        let lo = self.codes.partition_point(|&c| c < code);
        let hi = lo + self.codes[lo..].partition_point(|&c| c == code);
        &self.pos[lo..hi]
    }

    /// Each distinct code with its positions, in code order.
    pub(super) fn runs(&self) -> impl Iterator<Item = (u64, &[u32])> {
        let mut i = 0;
        std::iter::from_fn(move || {
            let code = *self.codes.get(i)?;
            let end = run_end(&self.codes, i);
            let run = (code, &self.pos[i..end]);
            i = end;
            Some(run)
        })
    }
}

/// End of the run of equal codes starting at `i`.
fn run_end(codes: &[u64], i: usize) -> usize {
    let mut end = i + 1;
    while end < codes.len() && codes[end] == codes[i] {
        end += 1;
    }
    end
}

/// Every (read position, backbone position) pair at which the two sequences share a
/// minimizer code: each shared code pairs each of its read positions with each of its
/// backbone positions. A merge of the two sorted code arrays, so no hashing.
pub(super) fn shared_anchors(read: &MinimizerMap, backbone: &MinimizerMap) -> Vec<(u32, u32)> {
    let (a, b) = (&read.codes, &backbone.codes);
    let (mut i, mut j) = (0, 0);
    let mut anchors = Vec::new();
    while i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            Ordering::Less => i += 1,
            Ordering::Greater => j += 1,
            Ordering::Equal => {
                let (i_end, j_end) = (run_end(a, i), run_end(b, j));
                for &p in &read.pos[i..i_end] {
                    for &q in &backbone.pos[j..j_end] {
                        anchors.push((p, q));
                    }
                }
                i = i_end;
                j = j_end;
            }
        }
    }
    anchors
}

/// 2-bit code and start position of every k-mer; k-mers containing a non-ACGT base are skipped.
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

/// Window-`w` minimizers of `seq`: each code with every position where it is the window
/// minimum. Most codes carry exactly one position; a code the window picks more than once
/// (a repeated minimizer) carries all of them.
pub(super) fn minimizer_map(seq: &[u8], k: usize, w: usize) -> MinimizerMap {
    let kms = encode_kmers(seq, k);
    let mut picks: Vec<(u64, u32)> = Vec::new();
    if kms.len() < w.max(1) {
        picks = kms;
    } else {
        let (mut last, mut lastpos) = (u64::MAX, u32::MAX);
        for win in kms.windows(w) {
            let m = *win.iter().min_by(|a, b| a.0.cmp(&b.0).then(a.1.cmp(&b.1))).unwrap();
            if m.0 != last || m.1 != lastpos {
                picks.push(m);
                last = m.0;
                lastpos = m.1;
            }
        }
    }
    picks.retain(|&(code, _)| !is_homopolymer(code, k));
    picks.sort_unstable();
    MinimizerMap {
        codes: picks.iter().map(|&(c, _)| c).collect(),
        pos: picks.iter().map(|&(_, p)| p).collect(),
    }
}

/// Whether every base of `code`'s k-mer is the same, as in a polyA tail. Such a k-mer is the
/// minimizer at every position of its run, so keeping all of them would pair each position of
/// one run with each of the other's: a ladder of equally valid anchors that the chain follows
/// to an arbitrary offset, putting the terminal anchor inside the tail instead of at the last
/// real sequence.
fn is_homopolymer(code: u64, k: usize) -> bool {
    let base = code & 3;
    (0..k).all(|i| (code >> (2 * i)) & 3 == base)
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
pub(super) fn chain_anchors(anchors: &mut [(u32, u32)]) -> Vec<(u32, u32)> {
    anchors.sort_unstable_by_key(|&(i, j)| (i, Reverse(j)));
    let chain = greedy_chain(anchors);
    if chain.len() == anchors.len() {
        return chain;
    }
    lis_chain(anchors)
}

/// Keep each anchor whose backbone position is above the last one kept. Anchors that share a
/// read position are sorted backbone-descending, so at most the first one seen can pass.
fn greedy_chain(anchors: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let mut chain: Vec<(u32, u32)> = Vec::with_capacity(anchors.len());
    for &(i, j) in anchors {
        match chain.last() {
            Some(&(_, lj)) if j <= lj => continue,
            _ => chain.push((i, j)),
        }
    }
    chain
}

/// Index buffers for [`lis_chain`], reused across calls on each worker thread because
/// `fits_interval` chains millions of anchor sets.
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
fn lis_chain(anchors: &[(u32, u32)]) -> Vec<(u32, u32)> {
    if anchors.is_empty() {
        return Vec::new();
    }
    SCRATCH.with(|s| {
        let s = &mut *s.borrow_mut();
        s.tails.clear();
        s.prev.clear();
        s.prev.resize(anchors.len(), NO_PREV);

        for (idx, &(_, j)) in anchors.iter().enumerate() {
            // Strictly increasing: the first tail whose backbone position is >= j.
            let pos = s.tails.partition_point(|&t| anchors[t as usize].1 < j);
            if pos > 0 {
                s.prev[idx] = s.tails[pos - 1];
            }
            if pos == s.tails.len() {
                s.tails.push(idx as u32);
            } else {
                s.tails[pos] = idx as u32;
            }
        }
        let mut out = Vec::with_capacity(s.tails.len());
        let mut cur = *s.tails.last().unwrap();
        loop {
            out.push(anchors[cur as usize]);
            if s.prev[cur as usize] == NO_PREV {
                break;
            }
            cur = s.prev[cur as usize];
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
        // "AAAC" (code 1) is the window minimum at both position 0 and position 10. The old
        // unique-only filter dropped such a code entirely; every occurrence is kept now.
        let map = minimizer_map(b"AAACTTTTTTAAACTTTTTT", 4, 2);
        assert_eq!(map.get(1), &[0u32, 10][..]);
    }

    #[test]
    fn shared_anchors_pairs_every_occurrence_of_each_shared_code() {
        // Code 1 at read 0 and 10 and backbone 5; code 7 at read 3 and backbone 8 and 20;
        // code 9 only in the read and code 4 only in the backbone.
        let read = MinimizerMap {
            codes: vec![1, 1, 7, 9].into(),
            pos: vec![0, 10, 3, 6].into(),
        };
        let backbone = MinimizerMap {
            codes: vec![1, 4, 7, 7].into(),
            pos: vec![5, 2, 8, 20].into(),
        };
        let mut anchors = shared_anchors(&read, &backbone);
        anchors.sort_unstable();
        assert_eq!(anchors, vec![(0, 5), (3, 8), (3, 20), (10, 5)]);
        assert_eq!(read.runs().map(|(c, p)| (c, p.len())).collect::<Vec<_>>(), vec![(1, 2), (7, 1), (9, 1)]);
    }

    #[test]
    fn minimizer_map_excludes_homopolymers() {
        // Every k-mer of a polyA run is the same code, one per position: a ladder of equally
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
        // A minimizer repeated in the read (i=30 vs two backbone hits) AND one repeated in
        // the backbone (j=100 vs two read hits) in the same anchor set.
        let mut anchors = vec![(10, 100), (30, 110), (30, 90), (50, 100)];
        let chain = chain_anchors(&mut anchors);
        assert_eq!(chain.len(), 2, "no length-3 collinear chain exists among these anchors");
        assert!(
            chain.windows(2).all(|w| w[0].0 < w[1].0 && w[0].1 < w[1].1),
            "chain must strictly increase in both coordinates: {chain:?}"
        );
    }
}
