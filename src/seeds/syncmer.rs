// ─────────────────────────────────────────────────────────────────────────────
// DUPLICATED FROM synpact/src/syncmer.rs (open-syncmer selection only).
// ─────────────────────────────────────────────────────────────────────────────
use super::hash::*;
use std::collections::VecDeque;

/// Lightweight syncmer record: position + the full-k-mer base-4 value.
#[derive(Clone)]
pub struct SyncmerLight {
    pub pos:   u32,
    pub value: u64,
}

/// Open-syncmer selection: a k-mer is a syncmer iff its minimum s-mer hash sits
/// at offset t. Linear-time via a monotone deque over the (k-s+1)-wide window.
pub fn select_syncmers_light(seq: &[u8], k: usize, s: usize, t: usize) -> Vec<SyncmerLight> {
    if seq.len() < k { return vec![]; }
    let mut hashes = match DnaHashFwd::new(seq, s) {
        Some(it) => it,
        None => return vec![],
    };
    let w = k - s + 1;
    let n_pos = seq.len() - k + 1;
    let mut out = Vec::with_capacity(n_pos / w + 1);
    let mut dq: VecDeque<(usize, u64)> = VecDeque::new();
    let mut r = 0usize;
    for i in 0..n_pos {
        let right = i + w - 1;
        while r <= right {
            let h = hashes.next().unwrap();
            while let Some(&(_, hb)) = dq.back() {
                if hb > h { dq.pop_back(); } else { break; }
            }
            dq.push_back((r, h));
            r += 1;
        }
        while let Some(&(f, _)) = dq.front() {
            if f < i { dq.pop_front(); } else { break; }
        }
        let min_j = dq.front().unwrap().0 - i;
        if min_j == t {
            let kmer = &seq[i..i + k];
            if kmer.iter().any(|&b| matches!(b.to_ascii_uppercase(), b'N')) { continue; }
            out.push(SyncmerLight { pos: i as u32, value: atom_value(kmer) });
        }
    }
    out
}
