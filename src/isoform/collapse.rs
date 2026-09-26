//! Step 5: merge isoforms that turn out to be the same transcript.
//!
//! Resolution can split one transcript into several isoforms: reads whose errors break
//! their anchors fail step 1 and gather in small satellite groups, and a read rejected by
//! one backbone can found a duplicate group. Consensuses are error-corrected, so they can be
//! compared more strictly than reads were.
//!
//! Isoforms are visited largest first, and each merges into the first larger one kept so far
//! that either
//!   * has a near-identical consensus (`collapse_ident`), whatever the read counts, or
//!   * contains it by the step-1 test with `max_gap = collapse_gap` and `min_cov >= 0.95`,
//!     if it has at most `collapse_ratio` of that isoform's reads. Comparable read counts
//!     suggest real isoforms with different ends rather than a satellite.
//!
//! A merged isoform's reads join the kept isoform, which keeps its own consensus.

use std::cmp::Reverse;
use std::collections::HashMap;

use super::anchors::{minimizer_map, MinimizerMap};
use super::fit::{fits_interval, strip_homopolymer_ends};
use super::{ba, Cfg, Isoform};

pub(super) fn collapse(mut isos: Vec<Isoform>, cfg: &Cfg) -> Vec<Isoform> {
    if isos.len() < 2 {
        return isos;
    }
    isos.sort_by_key(|iso| Reverse(iso.members.len()));
    let maps: Vec<MinimizerMap> =
        isos.iter().map(|iso| minimizer_map(&iso.consensus, cfg.k, cfg.w)).collect();
    let strict = Cfg { min_cov: cfg.min_cov.max(0.95), max_gap: cfg.collapse_gap, ..cfg.clone() };

    let n = isos.len();
    let mut target: Vec<Option<usize>> = vec![None; n]; // merged isoform -> kept isoform
    let mut kept: Vec<usize> = Vec::new();
    for i in 0..n {
        let into = kept.iter().copied().find(|&j| {
            should_merge(&isos[i], &isos[j], &maps[i], &maps[j], cfg, &strict)
        });
        match into {
            Some(j) => target[i] = Some(j),
            None => kept.push(i),
        }
    }

    // Kept isoforms in visiting order, each absorbing the reads merged into it.
    let mut pos: HashMap<usize, usize> = HashMap::new();
    let mut out: Vec<Isoform> = Vec::with_capacity(kept.len());
    for (p, &k) in kept.iter().enumerate() {
        pos.insert(k, p);
        out.push(Isoform {
            members: std::mem::take(&mut isos[k].members),
            consensus: std::mem::take(&mut isos[k].consensus),
        });
    }
    for i in 0..n {
        if let Some(j) = target[i] {
            let extra = std::mem::take(&mut isos[i].members);
            out[pos[&j]].members.extend(extra);
        }
    }
    out
}

/// Whether `small` (with no more reads than `big`) should merge into `big`. `strict` is
/// `cfg` tightened for comparing consensuses.
fn should_merge(
    small: &Isoform,
    big: &Isoform,
    small_map: &MinimizerMap,
    big_map: &MinimizerMap,
    cfg: &Cfg,
    strict: &Cfg,
) -> bool {
    if cfg.collapse_ident > 0.0 && same_sequence(small, big, cfg) {
        return true;
    }
    if small.members.len() as f64 > cfg.collapse_ratio * big.members.len() as f64 {
        return false;
    }
    // Containment, testing the shorter consensus against the longer.
    if small.consensus.len() <= big.consensus.len() {
        fits_interval(small_map, big_map, &small.consensus, &big.consensus, strict).is_some()
    } else {
        fits_interval(big_map, small_map, &big.consensus, &small.consensus, strict).is_some()
    }
}

/// Near-identical consensuses, ignoring homopolymer ends: identity at least `collapse_ident`,
/// with the longest indel and the length difference both within `collapse_gap`.
fn same_sequence(a: &Isoform, b: &Isoform, cfg: &Cfg) -> bool {
    let a = strip_homopolymer_ends(&a.consensus);
    let b = strip_homopolymer_ends(&b.consensus);
    let st = ba::gap_runs(a, b);
    st.ident() >= cfg.collapse_ident
        && st.max_indel <= cfg.collapse_gap as usize
        && (a.len() as i64 - b.len() as i64).unsigned_abs() as u32 <= cfg.collapse_gap
}
