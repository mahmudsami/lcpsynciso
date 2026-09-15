// ─────────────────────────────────────────────────────────────────────────────
// DUPLICATED FROM synpact/src/lcp.rs (parsing core only; no config/profiling).
// Adds `block_hashes_per_level` — the analog of synpact's `extract_all_levels`
// but it keeps *every* level's block hashes (no --min-level filtering) and drops
// positions/mass, since the histogram only needs the hashes.
// ─────────────────────────────────────────────────────────────────────────────
use super::hash::*;
use super::syncmer::*;

/// One LCP block: half-open range of unit indices `start..end`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Block {
    pub start: u32,
    pub end:   u32, // exclusive
}

impl Block {
    #[inline] pub fn range(&self) -> std::ops::Range<usize> {
        self.start as usize .. self.end as usize
    }
    #[inline] pub fn first(&self) -> usize { self.start as usize }
}

/// Cole–Vishkin deterministic coin-tossing → proper 3-colouring over {0,1,2}.
fn dct_three_coloring(vals: &[u64]) -> Vec<u8> {
    let m = vals.len();
    debug_assert!(m >= 2);
    let mut coin: Vec<u64> = vals.to_vec();

    while coin.iter().copied().max().unwrap_or(0) >= 6 {
        let mut next = vec![0u64; m];
        for i in 0..m - 1 {
            let p = (coin[i] ^ coin[i + 1]).trailing_zeros() as u64;
            next[i] = 2 * p + ((coin[i] >> p) & 1);
        }
        next[m - 1] = if next[m - 2] == 0 { 1 } else { 0 };
        coin = next;
    }

    let mut col: Vec<u8> = coin.iter().map(|&c| c as u8).collect();
    for v in [3u8, 4, 5] {
        for i in 0..m {
            if col[i] != v { continue; }
            let l = if i > 0     { col[i - 1] } else { u8::MAX };
            let r = if i + 1 < m { col[i + 1] } else { u8::MAX };
            col[i] = (0u8..=2).find(|&c| c != l && c != r).unwrap();
        }
    }
    col
}

/// Split a strictly-monotone run into ≤3-unit blocks, deterministically.
pub fn split_monotone_run(run: &[usize], values: &[u64]) -> Vec<(u32, u32)> {
    let m = run.len();
    if m <= 3 { return vec![(run[0] as u32, run[m - 1] as u32 + 1)]; }

    let vals: Vec<u64> = run.iter().map(|&i| values[i]).collect();
    let col: Vec<u64> = dct_three_coloring(&vals).iter().map(|&c| c as u64).collect();

    let triplet = |i: usize| -> (u32, u32) {
        let lo = i.saturating_sub(1);
        let hi = (i + 1).min(m - 1);
        (run[lo] as u32, run[hi] as u32 + 1)
    };

    let mut out: Vec<(u32, u32)> = Vec::new();
    for i in 0..m {
        if is_local_min(&col, i) { out.push(triplet(i)); }
    }
    for i in 0..m {
        if is_local_max(&col, i) {
            let l_min = i > 0     && is_local_min(&col, i - 1);
            let r_min = i + 1 < m && is_local_min(&col, i + 1);
            if !l_min && !r_min { out.push(triplet(i)); }
        }
    }
    out
}

#[inline]
pub fn is_local_min(v: &[u64], i: usize) -> bool {
    let n = v.len();
    (i == 0 || v[i] < v[i - 1]) && (i == n - 1 || v[i] < v[i + 1])
}

#[inline]
pub fn is_local_max(v: &[u64], i: usize) -> bool {
    let n = v.len();
    (i == 0 || v[i] > v[i - 1]) && (i == n - 1 || v[i] > v[i + 1])
}

/// Parse a stream of unit values into non-overlapping (boundary-sharing) blocks.
pub fn locally_consistent_parsing(values: &[u64]) -> Vec<Block> {
    let n = values.len();
    if n == 0 { return vec![]; }
    if n == 1 { return vec![Block { start: 0, end: 1 }]; }

    let is_min: Vec<bool> = (0..n).map(|i| is_local_min(values, i)).collect();
    let is_max: Vec<bool> = (0..n).map(|i| is_local_max(values, i)).collect();

    let mut assigned = vec![false; n];
    let mut blocks: Vec<Block> = Vec::new();

    macro_rules! commit_range {
        ($lo:expr, $hi:expr) => {{
            let lo = $lo;
            let hi = ($hi).min(n - 1);
            if (lo..=hi).any(|i| !assigned[i]) {
                for i in lo..=hi { assigned[i] = true; }
                blocks.push(Block { start: lo as u32, end: hi as u32 + 1 });
            }
        }};
    }

    // Rule 1 — local-minimum triplet.
    for i in 0..n {
        if is_min[i] { commit_range!(i.saturating_sub(1), i + 1); }
    }
    // Rule 2 — local-maximum triplet (no adjacent local minimum).
    for i in 0..n {
        if is_max[i] {
            let l_is_min = i > 0     && is_min[i - 1];
            let r_is_min = i + 1 < n && is_min[i + 1];
            if !l_is_min && !r_is_min {
                commit_range!(i.saturating_sub(1), i + 1);
            }
        }
    }
    // Rule 3 — repetition run ≥ 2 plus neighbours.
    {
        let mut i = 0;
        while i < n {
            let mut j = i + 1;
            while j < n && values[j] == values[i] { j += 1; }
            if j - i >= 2 { commit_range!(i.saturating_sub(1), j); }
            i = j;
        }
    }
    // Rule 4 — full monotone run, split via DCT.
    for a in 0..n.saturating_sub(2) {
        if !assigned[a] { continue; }
        for &sign in &[1i128, -1i128] {
            let mut b = a + 1;
            while b < n {
                let d = values[b] as i128 - values[b - 1] as i128;
                if d * sign > 0 { b += 1; } else { break; }
            }
            b -= 1;
            if b > a + 1 && assigned[b] && (a + 1..b).any(|j| !assigned[j]) {
                let run: Vec<usize> = (a..=b).collect();
                for (s, e) in split_monotone_run(&run, values) {
                    commit_range!(s as usize, e as usize - 1);
                }
            }
        }
    }

    let bad: Vec<usize> = (0..n).filter(|&i| !assigned[i]).collect();
    if !bad.is_empty() {
        panic!(
            "LCP bug: {} position(s) unassigned (first few: {:?})\nvalues={:?}",
            bad.len(), &bad[..bad.len().min(5)], values
        );
    }

    blocks.sort_by_key(|b| b.start);
    blocks
}

/// For one read sequence, return the block hashes at every level 1..=max_levels.
/// `out[0]` = L1 hashes, `out[1]` = L2, … Mirrors synpact's level recursion
/// (`extract_all_levels`) but with no min-level filtering and no positions.
pub fn block_hashes_per_level(
    seq: &[u8], k: usize, s: usize, t: usize, max_levels: usize,
) -> Vec<Vec<u64>> {
    let syncmers = select_syncmers_light(seq, k, s, t);
    if syncmers.is_empty() { return vec![]; }

    let smer_vals: Vec<u64> = syncmers.iter().map(|sm| sm.value).collect();
    let l1_raw = locally_consistent_parsing(&smer_vals);
    if l1_raw.is_empty() { return vec![]; }

    let mut out: Vec<Vec<u64>> = Vec::with_capacity(max_levels);

    // L1 block hashes.
    let mut cur_hashes: Vec<u64> = Vec::with_capacity(l1_raw.len());
    for blk in &l1_raw {
        cur_hashes.push(block_hash_for_level(&smer_vals[blk.range()], 0));
    }
    out.push(cur_hashes.clone());

    // L2 … Lmax: parse the level below and re-hash.
    for level_1idx in 2..=max_levels {
        if cur_hashes.len() < 2 { break; }
        let next_raw = locally_consistent_parsing(&cur_hashes);
        if next_raw.is_empty() { break; }

        let mut next_hashes = Vec::with_capacity(next_raw.len());
        for blk in &next_raw {
            next_hashes.push(block_hash_for_level(&cur_hashes[blk.range()], level_1idx - 1));
        }
        out.push(next_hashes.clone());
        cur_hashes = next_hashes;
    }

    out
}
