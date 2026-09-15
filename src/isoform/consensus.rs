//! Step 4: consensus per isoform, from anchor-windowed modal segments.
//!
//! The reads of one isoform are near-identical and collinear (that is exactly
//! why they were grouped). So instead of a multiple-sequence alignment / POA we
//! exploit the shared unique minimizer anchors we already know how to compute:
//!
//!   * backbone      = the longest read;
//!   * WINDOW BOUNDARIES are the backbone minimizers that are ALSO shared by a
//!     majority of the reads — an anchor a single read's error can't destroy.
//!     (Backbone-unique minimizers are skipped: they cluster exactly at the
//!     backbone's own errors, and anchoring there would leave that error with
//!     no other read able to vote it out — the whole point is to overrule it.)
//!   * a read contributes to a window iff it shares the minimizers at BOTH of
//!     the window's flanking backbone positions (a collinear anchor pair) — its
//!     substring between those two anchors is that window's candidate;
//!   * the window's consensus is the MOST FREQUENT exact candidate substring, so
//!     a backbone error sits inside a majority-anchored window and is outvoted.
//!
//! Because the flanking k-mers are exact matches in every contributing read, the
//! segments are perfectly registered and concatenate seamlessly — no inner
//! alignment at all. Cost is ~O(total read length): one minimizer pass per read
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

use super::Read;

/// Refine an isoform's consensus from its member reads (indices into `reads`),
/// reusing minimizer maps already computed for the whole cluster (`maps[ri]`).
/// Falls back to the backbone sequence for isoforms too small or too
/// minimizer-poor to vote on.
pub fn refine_consensus_with_maps(
    reads: &[Read],
    members: &[usize],
    maps: &[HashMap<u64, u32>],
) -> Vec<u8> {
    // Backbone = longest member.
    let bb = *members.iter().max_by_key(|&&i| reads[i].seq.len()).unwrap();
    let backbone = &reads[bb].seq;
    if members.len() < 2 {
        return backbone.clone();
    }

    // How many members carry each minimizer code.
    let mut freq: HashMap<u64, u32> = HashMap::new();
    for &ri in members {
        for &c in maps[ri].keys() {
            *freq.entry(c).or_insert(0) += 1;
        }
    }
    let thresh = (members.len() / 2 + 1) as u32; // strict majority

    // Breakpoints = backbone minimizers shared by a majority, ordered by
    // position. Positions are k-mer starts, all distinct (one code per
    // position), so the sort is strictly increasing.
    let mut breaks: Vec<(u32, u64)> = maps[bb]
        .iter()
        .filter(|(c, _)| freq[*c] >= thresh)
        .map(|(&c, &p)| (p, c))
        .collect();
    breaks.sort_unstable();
    if breaks.len() < 2 {
        return backbone.clone();
    }
    let m = breaks.len();
    let n_windows = m - 1;

    // Per-window vote table: candidate segment bytes -> count.
    let mut votes: Vec<HashMap<Vec<u8>, u32>> = vec![HashMap::new(); n_windows];

    for &ri in members {
        let seq = &reads[ri].seq;
        let rmap = &maps[ri];
        // Read position at each backbone breakpoint, kept as a strictly
        // increasing (collinear) chain; non-matching or out-of-order -> None.
        let mut last = -1i64;
        let mut rpos: Vec<Option<u32>> = Vec::with_capacity(m);
        for &(_p, c) in &breaks {
            match rmap.get(&c) {
                Some(&rp) if (rp as i64) > last => {
                    rpos.push(Some(rp));
                    last = rp as i64;
                }
                _ => rpos.push(None),
            }
        }
        for wi in 0..n_windows {
            if let (Some(a), Some(b)) = (rpos[wi], rpos[wi + 1]) {
                let seg = seq[a as usize..b as usize].to_vec();
                *votes[wi].entry(seg).or_insert(0) += 1;
            }
        }
    }

    // Assemble: unanchored leading end, modal windows, unanchored trailing end.
    let mut out: Vec<u8> = Vec::with_capacity(backbone.len());
    out.extend_from_slice(&backbone[..breaks[0].0 as usize]);
    for wi in 0..n_windows {
        let default = &backbone[breaks[wi].0 as usize..breaks[wi + 1].0 as usize];
        out.extend_from_slice(&modal_segment(&votes[wi], default));
    }
    out.extend_from_slice(&backbone[breaks[m - 1].0 as usize..]);
    out
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
