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
//! A merged isoform's reads join the kept isoform, which keeps its own consensus, unless the
//! merged one is a more complete copy of it — reaching further 5' than a 5'-truncated
//! majority, or on to the real 3' end past reads that stop early (internal priming): then the
//! kept isoform's consensus is rebuilt from all its reads, so the full-length sequence survives.
//! Containment treats an indel of more than `collapse_gap` bp as a difference, as the
//! near-identity test does.

use std::cmp::Reverse;
use std::collections::HashMap;

use super::anchors::{minimizer_map, MinimizerMap};
use super::fit::{fits_interval, strip_homopolymer_ends};
use super::{ba, consensus, Cfg, Isoform, Read};

/// Resolve duplicates among `isos`, whose members index `reads` (with minimizer maps
/// `read_maps`, used to rebuild a consensus).
pub(super) fn collapse(
    mut isos: Vec<Isoform>,
    reads: &[Read],
    read_maps: &[MinimizerMap],
    cfg: &Cfg,
) -> Vec<Isoform> {
    if isos.len() < 2 {
        return isos;
    }
    isos.sort_by_key(|iso| Reverse(iso.members.len()));
    let maps: Vec<MinimizerMap> =
        isos.iter().map(|iso| minimizer_map(&iso.consensus, cfg.k, cfg.w)).collect();
    let strict = Cfg {
        min_cov: cfg.min_cov.max(0.95),
        max_gap: cfg.collapse_gap,
        max_indel_run: cfg.collapse_gap + 1,
        ..cfg.clone()
    };

    let n = isos.len();
    let mut target: Vec<Option<usize>> = vec![None; n]; // merged isoform -> kept isoform
    let mut rebuild = vec![false; n]; // kept isoform absorbed a more complete copy of itself
    let mut kept: Vec<usize> = Vec::new();
    for i in 0..n {
        let mut into = None;
        for &j in &kept {
            match should_merge(&isos[i], &isos[j], &maps[i], &maps[j], cfg, &strict) {
                Merge::No => continue,
                Merge::Yes => {}
                Merge::Extends => rebuild[j] = true,
            }
            into = Some(j);
            break;
        }
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
            alt_start: isos[k].alt_start,
        });
    }
    for i in 0..n {
        if let Some(j) = target[i] {
            let extra = std::mem::take(&mut isos[i].members);
            out[pos[&j]].members.extend(extra);
        }
    }
    for (p, &k) in kept.iter().enumerate() {
        if rebuild[k] {
            let iso = &mut out[p];
            iso.consensus = if cfg.consensus {
                consensus::refine_consensus_with_maps(reads, &iso.members, &iso.members, read_maps)
            } else {
                let longest = *iso.members.iter().max_by_key(|&&r| reads[r].seq.len()).unwrap();
                reads[longest].seq.clone()
            };
        }
    }
    out
}

/// How `small` relates to a larger isoform it may merge into.
enum Merge {
    No,
    Yes,
    /// Merge, and `small` is a more complete copy of `big`: rebuild the consensus.
    Extends,
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
) -> Merge {
    if cfg.collapse_ident > 0.0 && same_sequence(small, big, cfg) {
        return Merge::Yes;
    }
    if small.members.len() as f64 > cfg.collapse_ratio * big.members.len() as f64 {
        return Merge::No;
    }
    // Containment, testing the shorter consensus against the longer. An end reaching more
    // than `collapse_gap` bp past the other's is a real extension, not a consensus wobble.
    // An isoform at a validated downstream start is a transcript of its own, so a pair that
    // differs at the 5' end stays apart when either is one.
    let gap = cfg.collapse_gap as i32;
    let alt = small.alt_start || big.alt_start;
    if small.consensus.len() <= big.consensus.len() {
        // `small` placed on `big`; it may overhang `big`'s 5' end (a negative start).
        match fits_interval(small_map, big_map, &small.consensus, &big.consensus, strict) {
            None => Merge::No,
            Some((s, e)) if s < -gap => {
                // Same 3' end: `small` is `big` with its 5' end intact. Otherwise the two
                // differ at both ends and stay apart.
                if !alt && e >= big.consensus.len() as i32 - cfg.boundary_tol as i32 {
                    Merge::Extends
                } else {
                    Merge::No
                }
            }
            Some(_) => Merge::Yes,
        }
    } else {
        // `big` placed on `small`: `small` extends it at the 5' end, the 3' end, or both.
        match fits_interval(big_map, small_map, &big.consensus, &small.consensus, strict) {
            None => Merge::No,
            Some((s, _)) if s > gap && alt => Merge::No,
            Some((s, e)) if s > gap || e < small.consensus.len() as i32 - gap => Merge::Extends,
            Some(_) => Merge::Yes,
        }
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

#[cfg(test)]
mod tests {
    use super::*;

    fn random_seq(n: usize, seed: u64) -> Vec<u8> {
        let mut x = seed;
        (0..n)
            .map(|_| {
                x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                b"ACGT"[(x >> 33) as usize % 4]
            })
            .collect()
    }

    #[test]
    fn a_5_prime_complete_minority_gives_the_merged_isoform_its_full_sequence() {
        // Ten reads missing the first 500 bp form the larger isoform; three full-length reads
        // form a smaller one. They merge, and the kept isoform must not stay truncated.
        let truth = random_seq(1500, 5);
        let mut seqs: Vec<Vec<u8>> = (0..10).map(|_| truth[500..].to_vec()).collect();
        seqs.extend((0..3).map(|_| truth.clone()));
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| Read { name: i.to_string(), seq: s.clone() })
            .collect();
        let cfg = Cfg::default();
        let maps: Vec<MinimizerMap> =
            reads.iter().map(|r| minimizer_map(&r.seq, cfg.k, cfg.w)).collect();
        let isos = vec![
            Isoform {
                members: (0..10).collect(),
                consensus: truth[500..].to_vec(),
                alt_start: false,
            },
            Isoform {
                members: (10..13).collect(),
                consensus: truth.clone(),
                alt_start: false,
            },
        ];
        let out = collapse(isos, &reads, &maps, &cfg);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].members.len(), 13);
        assert_eq!(out[0].consensus, truth);
    }

    #[test]
    fn an_isoform_at_a_validated_start_is_not_folded_into_the_full_length_one() {
        let truth = random_seq(1500, 5);
        let mut seqs: Vec<Vec<u8>> = (0..10).map(|_| truth[500..].to_vec()).collect();
        seqs.extend((0..3).map(|_| truth.clone()));
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| Read { name: i.to_string(), seq: s.clone() })
            .collect();
        let cfg = Cfg::default();
        let maps: Vec<MinimizerMap> =
            reads.iter().map(|r| minimizer_map(&r.seq, cfg.k, cfg.w)).collect();
        let isos = vec![
            Isoform {
                members: (0..10).collect(),
                consensus: truth[500..].to_vec(),
                alt_start: true,
            },
            Isoform {
                members: (10..13).collect(),
                consensus: truth.clone(),
                alt_start: false,
            },
        ];
        assert_eq!(collapse(isos, &reads, &maps, &cfg).len(), 2);
    }

    #[test]
    fn reads_reaching_the_real_3_prime_end_give_the_merged_isoform_its_full_sequence() {
        // Ten reads stop 600 bp early (internal priming) and form the larger isoform; three
        // reach the real 3' end. The merged isoform must keep the full-length sequence.
        let truth = random_seq(1500, 9);
        let mut seqs: Vec<Vec<u8>> = (0..10).map(|_| truth[..900].to_vec()).collect();
        seqs.extend((0..3).map(|_| truth.clone()));
        let reads: Vec<Read> = seqs
            .iter()
            .enumerate()
            .map(|(i, s)| Read { name: i.to_string(), seq: s.clone() })
            .collect();
        let cfg = Cfg::default();
        let maps: Vec<MinimizerMap> =
            reads.iter().map(|r| minimizer_map(&r.seq, cfg.k, cfg.w)).collect();
        let isos = vec![
            Isoform {
                members: (0..10).collect(),
                consensus: truth[..900].to_vec(),
                alt_start: false,
            },
            Isoform {
                members: (10..13).collect(),
                consensus: truth.clone(),
                alt_start: false,
            },
        ];
        let out = collapse(isos, &reads, &maps, &cfg);
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].consensus, truth);
    }
}
