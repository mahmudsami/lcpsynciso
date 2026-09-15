//! Minimizer anchors between two sequences, and the collinear chain through them.
//!
//! An anchor is a minimizer that occurs once in each sequence, pairing a read position with
//! a backbone position. The chain is the largest set of anchors whose positions increase
//! together in both sequences.

use std::cell::RefCell;
use std::collections::HashMap;

// ── Minimizers ───────────────────────────────────────────────────────────────

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

/// Window-`w` minimizers of `seq`, as code -> position. A code that is the minimizer at more
/// than one position is left out, so every anchor built from the map is unambiguous.
pub(super) fn minimizer_map(seq: &[u8], k: usize, w: usize) -> HashMap<u64, u32> {
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
    // Position per code, or -1 once a code has been seen at a second position.
    let mut seen: HashMap<u64, i64> = HashMap::with_capacity(picks.len());
    for (code, pos) in picks {
        seen.entry(code).and_modify(|v| *v = -1).or_insert(pos as i64);
    }
    seen.into_iter().filter(|&(_, v)| v >= 0).map(|(c, v)| (c, v as u32)).collect()
}

// ── Chaining ─────────────────────────────────────────────────────────────────

/// The longest chain of `anchors` (sorted by read position) whose backbone positions strictly
/// increase.
///
/// Maximising the number of anchors keeps a mispaired anchor out of the chain, because
/// keeping it would block every anchor it sits in front of. Most comparisons have no
/// mispaired anchor, so the greedy scan runs first: when it keeps every anchor, the full set
/// is already the longest chain.
pub(super) fn chain_anchors(anchors: &[(u32, u32)]) -> Vec<(u32, u32)> {
    let chain = greedy_chain(anchors);
    if chain.len() == anchors.len() {
        return chain;
    }
    lis_chain(anchors)
}

/// Keep each anchor whose backbone position is above the last one kept.
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
/// sorting with predecessor links).
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
