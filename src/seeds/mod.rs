//! LCP-syncmer seed machinery: open-syncmer selection, locally-consistent
//! parsing, and per-level block hashing. Duplicated from synpact (kept in sync
//! so seed/block hashes match exactly) — do not edit synpact's own tree.

pub mod hash;
pub mod lcp;
pub mod syncmer;

pub use lcp::block_hashes_per_level;

/// Forward-strand `(seed hash, level)` list for a read, deduplicated by hash.
///
/// Level is 1-based (`L1..=levels`). Block hashes are domain-separated per
/// level, so a hash belongs to exactly one level. FLNC reads are oriented
/// 5'->3', so same-transcript reads share forward-strand seeds — forward
/// strand only, no reverse complement.
pub fn read_seed_levels(
    seq: &[u8],
    k: usize,
    s: usize,
    t: usize,
    levels: usize,
) -> Vec<(u64, u8)> {
    let per = block_hashes_per_level(seq, k, s, t, levels);
    let mut v: Vec<(u64, u8)> = Vec::new();
    for (i, lvl) in per.into_iter().enumerate() {
        let lv = (i + 1) as u8;
        for h in lvl {
            v.push((h, lv));
        }
    }
    v.sort_unstable_by_key(|x| x.0);
    v.dedup_by_key(|x| x.0);
    v
}
