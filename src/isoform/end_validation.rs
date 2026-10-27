//! Step 5: validate each isoform's 3' end by its polyadenylation signal.
//!
//! Cleavage and polyadenylation happen 10-30 nt downstream of a polyadenylation signal (PAS),
//! the hexamer AAUAAA or a close variant, so a real 3' end carries one just upstream. An end
//! where the oligo-dT primer annealed to an A-rich stretch inside the RNA (internal priming)
//! lies wherever such a stretch is, and carries a PAS there only by chance. On UHRR PacBio HiFi
//! (lcpsynciso min-iso 3 predictions, labelled by SQANTI3: a PolyASite site, or > 60% A
//! downstream and no site), a PAS starting 10-34 bp before the end marks 91.5% of real ends and
//! 15.9% of internally primed ones: the bins where real ends carry it at least twice as often
//! as internally primed ones, which carry it at the chance rate. It needs no genome, only the
//! consensus and the hexamers.
//!
//! An isoform whose end has no PAS merges into the isoform of its cluster that contains it
//! (the strict fit of step 4) and runs on to a validated end, the nearest such end:
//!   * if it starts at a validated downstream start that differs from the container's, it
//!     stays an isoform of its own, running from that start to the container's end;
//!   * otherwise its reads join the container, which takes its 5' part if that reaches
//!     further 5' (the container's own 5' end being then a truncated copy).
//! One that nothing contains keeps its end unvalidated and goes to the secondary output files.

use super::anchors::{minimizer_map, MinimizerMap};
use super::collapse::strict_options;
use super::fit::{place_on_backbone, strip_homopolymer_ends};
use super::{stats, EndValidation, Isoform, IsoformOptions};

/// Polyadenylation-signal hexamers: AAUAAA and its common variants, the list SQANTI3 tests.
const PAS_HEXAMERS: [&[u8; 6]; 12] = [
    b"AATAAA", b"ATTAAA", b"AGTAAA", b"TATAAA", b"CATAAA", b"GATAAA", b"AATATA", b"AATACA",
    b"AATAGA", b"ACTAAA", b"AAGAAA", b"AATGAA",
];
/// A PAS starting this many bases before the 3' end (nearest, farthest) marks a real end.
const PAS_WINDOW: (usize, usize) = (10, 34);

/// Whether a PAS hexamer starts `PAS_WINDOW` bases before the 3' end of `consensus`, a
/// poly(A) tail (a homopolymer run of more than 10 bases, as ONT reads keep) not counted.
pub(super) fn has_pas(consensus: &[u8]) -> bool {
    let seq = strip_homopolymer_ends(consensus);
    let end = seq.len();
    let (near, far) = PAS_WINDOW;
    if end < near + PAS_HEXAMERS[0].len() {
        return false;
    }
    (end.saturating_sub(far)..=end - near).any(|i| {
        i + 6 <= end && PAS_HEXAMERS.iter().any(|hexamer| seq[i..i + 6].eq_ignore_ascii_case(&hexamer[..]))
    })
}

/// Validate the 3' ends of one cluster's isoforms (after merging), in their order.
pub(super) fn validate_ends(mut isoforms: Vec<Isoform>, options: &IsoformOptions) -> Vec<Isoform> {
    if options.end_validation == EndValidation::Off {
        return isoforms;
    }
    let n = isoforms.len();
    let signal: Vec<bool> = isoforms.iter().map(|isoform| has_pas(&isoform.consensus)).collect();
    let validated = signal.iter().filter(|&&s| s).count();
    stats::N_END_VALIDATED.fetch_add(validated as u64, std::sync::atomic::Ordering::Relaxed);
    if validated == n {
        return isoforms;
    }
    let maps: Vec<MinimizerMap> = isoforms
        .iter()
        .map(|isoform| minimizer_map(&isoform.consensus, options.minimizer_k, options.minimizer_w))
        .collect();
    let strict = strict_options(options);
    let slack = options.collapse_gap as i32;

    // Each isoform without a signal: the container with the nearest validated end, where on it
    // the isoform starts, and the container's sequence past the isoform's end.
    let mut container: Vec<Option<(usize, i32, Vec<u8>)>> = vec![None; n];
    for i in (0..n).filter(|&i| !signal[i]) {
        let mut best: Option<(usize, i32, i32, i32)> = None; // (container, start, end, extension)
        for j in (0..n).filter(|&j| j != i && signal[j]) {
            let Some((start, end)) = place_on_backbone(
                &maps[i],
                &maps[j],
                &isoforms[i].consensus,
                &isoforms[j].consensus,
                &strict,
            ) else {
                continue;
            };
            let extension = isoforms[j].consensus.len() as i32 - end;
            let nearer = |(b, _, _, e): (usize, i32, i32, i32)| {
                extension < e
                    || (extension == e && isoforms[j].members.len() > isoforms[b].members.len())
            };
            if extension > slack && best.map_or(true, nearer) {
                best = Some((j, start, end, extension));
            }
        }
        container[i] = best.map(|(j, start, end, _)| (j, start, isoforms[j].consensus[end as usize..].to_vec()));
    }

    let mut merged = vec![false; n];
    let mut took_5p = vec![false; n]; // a container's consensus already took a merged 5' part
    for i in (0..n).filter(|&i| !signal[i]) {
        let Some((j, start, tail)) = container[i].take() else {
            isoforms[i].end_validated = false;
            stats::inc(&stats::N_END_SECONDARY);
            continue;
        };
        if isoforms[i].alt_start && start.abs() > slack {
            // Its own validated start, and the container's validated end.
            isoforms[i].consensus.extend_from_slice(&tail);
            stats::inc(&stats::N_END_RESPLIT);
            continue;
        }
        let members = std::mem::take(&mut isoforms[i].members);
        if start < -slack && !took_5p[j] {
            let mut consensus = std::mem::take(&mut isoforms[i].consensus);
            consensus.extend_from_slice(&tail);
            isoforms[j].consensus = consensus;
            took_5p[j] = true;
        }
        isoforms[j].members.extend(members);
        merged[i] = true;
        stats::inc(&stats::N_END_MERGED);
    }
    isoforms.into_iter().zip(merged).filter(|&(_, gone)| !gone).map(|(isoform, _)| isoform).collect()
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

    /// `seq` with AATAAA starting `before` bases before its end.
    fn with_pas(mut seq: Vec<u8>, before: usize) -> Vec<u8> {
        let at = seq.len() - before;
        seq[at..at + 6].copy_from_slice(b"AATAAA");
        seq
    }

    /// A random sequence with no PAS hexamer anywhere in its last 60 bases.
    fn without_pas(n: usize, seed: u64) -> Vec<u8> {
        let mut seq = random_seq(n, seed);
        for b in seq[n - 60..].iter_mut() {
            if *b == b'A' {
                *b = b'C';
            }
        }
        seq
    }

    fn isoform(consensus: Vec<u8>, reads: std::ops::Range<usize>, alt_start: bool) -> Isoform {
        Isoform { members: reads.collect(), consensus, alt_start, end_validated: true }
    }

    /// A 1500 bp transcript ending with a PAS 22 bp before its end, with no PAS anywhere
    /// before position 1200 near the copies' ends (no A in 1140-1200).
    fn transcript() -> Vec<u8> {
        let mut seq = without_pas(1500, 7);
        for b in seq[1140..1200].iter_mut() {
            if *b == b'A' {
                *b = b'C';
            }
        }
        with_pas(seq, 22)
    }

    #[test]
    fn a_pas_counts_only_in_its_window_and_before_a_poly_a_tail() {
        let seq = without_pas(400, 3);
        assert!(!has_pas(&seq));
        assert!(has_pas(&with_pas(seq.clone(), 22)));
        assert!(!has_pas(&with_pas(seq.clone(), 50)));
        let mut tailed = with_pas(seq, 22);
        tailed.extend_from_slice(&[b'A'; 30]);
        assert!(has_pas(&tailed));
    }

    #[test]
    fn an_end_without_a_signal_merges_into_the_isoform_that_runs_on_to_a_validated_end() {
        // The full transcript ends with a PAS; a copy stopping 300 bp early has none.
        let full = transcript();
        let early = full[..1200].to_vec();
        assert!(has_pas(&full) && !has_pas(&early));
        let out = validate_ends(
            vec![isoform(full.clone(), 0..10, false), isoform(early, 10..14, false)],
            &IsoformOptions::default(),
        );
        assert_eq!(out.len(), 1);
        assert_eq!(out[0].members.len(), 14);
        assert!(out[0].end_validated);
        assert_eq!(out[0].consensus, full);
    }

    #[test]
    fn at_a_validated_start_it_stays_apart_and_takes_the_containers_end() {
        let full = transcript();
        let downstream = full[400..1200].to_vec();
        let n = downstream.len();
        assert!(!has_pas(&downstream));
        let out = validate_ends(
            vec![isoform(full.clone(), 0..10, false), isoform(downstream.clone(), 10..14, true)],
            &IsoformOptions::default(),
        );
        assert_eq!(out.len(), 2);
        let own = out.iter().find(|i| i.alt_start).unwrap();
        assert_eq!(own.members.len(), 4);
        assert!(own.end_validated);
        assert_eq!(&own.consensus[..n], &downstream[..]);
        assert_eq!(own.consensus.len(), 1100, "running on to the container's end");
    }

    #[test]
    fn an_end_without_a_signal_and_nothing_to_merge_into_is_kept_apart() {
        let out = validate_ends(
            vec![isoform(transcript(), 0..10, false), isoform(without_pas(900, 11), 10..14, false)],
            &IsoformOptions::default(),
        );
        assert_eq!(out.len(), 2);
        assert!(out[0].end_validated);
        assert!(!out[1].end_validated);
        let off = IsoformOptions { end_validation: EndValidation::Off, ..IsoformOptions::default() };
        assert!(validate_ends(vec![isoform(without_pas(900, 11), 0..4, false)], &off)[0].end_validated);
    }
}
