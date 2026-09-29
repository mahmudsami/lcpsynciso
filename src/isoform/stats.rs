//! Diagnostic counters for isoform detection, printed by `predict` after phase C.
//!
//! Counters are process-wide and summed over all worker threads, so times are CPU time in
//! nanoseconds, not wall time.

use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
use std::time::Instant;

// Time spent in each step of `detect_isoforms`.
pub(crate) static T_MAPS: AtomicU64 = AtomicU64::new(0);
pub(crate) static T_PARTITION: AtomicU64 = AtomicU64::new(0);
pub(crate) static T_CONS: AtomicU64 = AtomicU64::new(0);

// Inside `place_on_backbone`: time per stage, then event counts.
pub(crate) static T_ANCHOR: AtomicU64 = AtomicU64::new(0);
pub(crate) static T_GAP: AtomicU64 = AtomicU64::new(0);
pub(crate) static T_FLANK: AtomicU64 = AtomicU64::new(0);
/// Calls to `place_on_backbone`.
pub(crate) static N_PLACEMENT_CALLS: AtomicU64 = AtomicU64::new(0);
/// Interior gaps wide enough to be aligned.
pub(crate) static N_GAPS_ALIGNED: AtomicU64 = AtomicU64::new(0);
/// Traced interior-gap alignments.
pub(crate) static N_TRACE: AtomicU64 = AtomicU64::new(0);
/// Aligned interior gaps accepted as sequencing error.
pub(crate) static N_GAPS_ACCEPTED: AtomicU64 = AtomicU64::new(0);
/// Rejections, indexed by [`RejectReason`].
static N_REJ: [AtomicU64; 6] = [const { AtomicU64::new(0) }; 6];
/// Size histograms of aligned interior gaps and of flanks, bucketed by [`bucket`].
static H_GAP: [AtomicU64; 6] = [const { AtomicU64::new(0) }; 6];
static H_FLANK: [AtomicU64; 6] = [const { AtomicU64::new(0) }; 6];

// How `ends::split_by_ends` placed each read (reported as "split_by_ends routes").
/// Read fell in a dense (start, end) cell.
pub(crate) static N_DENSE: AtomicU64 = AtomicU64::new(0);
/// Read was folded into an isoform whose extent contains it.
pub(crate) static N_ENCLOSED: AtomicU64 = AtomicU64::new(0);
/// Read matched no isoform and was dropped.
pub(crate) static N_DROP: AtomicU64 = AtomicU64::new(0);
/// Folded read was wider than the extent it was folded into (expected to stay zero).
pub(crate) static N_ENCLOSED_WIDER: AtomicU64 = AtomicU64::new(0);

/// Why `place_on_backbone` rejected a read, in report order.
#[derive(Clone, Copy)]
pub(crate) enum RejectReason {
    FewAnchors,
    ShortChain,
    InteriorGap,
    FlankDiverged,
    BothFlanks,
    LowCoverage,
}

#[inline]
pub(crate) fn inc(counter: &AtomicU64) {
    counter.fetch_add(1, Relaxed);
}

#[inline]
pub(crate) fn add_elapsed(counter: &AtomicU64, since: Instant) {
    counter.fetch_add(since.elapsed().as_nanos() as u64, Relaxed);
}

#[inline]
pub(crate) fn count_rejection(reason: RejectReason) {
    inc(&N_REJ[reason as usize]);
}

#[inline]
pub(crate) fn record_gap_size(len: usize) {
    inc(&H_GAP[bucket(len)]);
}

#[inline]
pub(crate) fn record_flank_size(len: usize) {
    inc(&H_FLANK[bucket(len)]);
}

/// Histogram bucket: <32, 32-63, 64-127, 128-255, 256-511, >=512.
#[inline]
fn bucket(n: usize) -> usize {
    match n {
        0..=31 => 0,
        32..=63 => 1,
        64..=127 => 2,
        128..=255 => 3,
        256..=511 => 4,
        _ => 5,
    }
}

fn load(counters: &[AtomicU64; 6]) -> [u64; 6] {
    std::array::from_fn(|i| counters[i].load(Relaxed))
}

/// Nanoseconds spent so far on (minimizer maps, partition, consensus).
pub(crate) fn timings_ns() -> (u64, u64, u64) {
    (T_MAPS.load(Relaxed), T_PARTITION.load(Relaxed), T_CONS.load(Relaxed))
}

/// Print gap and flank size histograms, the `place_on_backbone` breakdown, how reads were placed
/// by end, and why reads were rejected.
pub(crate) fn print_report() {
    const BUCKETS: [&str; 6] = ["<32", "32-63", "64-127", "128-255", "256-511", ">=512"];
    const REJECT_LABELS: [&str; 6] =
        ["few-anchors", "short-chain", "INTERIOR-GAP", "flank-diverged", "both-flanks", "low-coverage"];
    const NANOS_PER_SEC: f64 = 1e9;

    let percentages = |h: [u64; 6]| {
        let total = h.iter().sum::<u64>().max(1);
        (0..6)
            .map(|i| format!("{}={:.1}%", BUCKETS[i], 100.0 * h[i] as f64 / total as f64))
            .collect::<Vec<_>>()
            .join("  ")
    };
    eprintln!("[predict] === gap sizes:   {} ===", percentages(load(&H_GAP)));
    eprintln!("[predict] === flank sizes: {} ===", percentages(load(&H_FLANK)));
    eprintln!(
        "[predict] === place_on_backbone breakdown: anchors+chain={:.0}s  interior-gap-DP={:.0}s  flank-DP={:.0}s  | calls={}  gap-DPs={}  traced={}  rescued-as-error={} ===",
        T_ANCHOR.load(Relaxed) as f64 / NANOS_PER_SEC,
        T_GAP.load(Relaxed) as f64 / NANOS_PER_SEC,
        T_FLANK.load(Relaxed) as f64 / NANOS_PER_SEC,
        N_PLACEMENT_CALLS.load(Relaxed),
        N_GAPS_ALIGNED.load(Relaxed),
        N_TRACE.load(Relaxed),
        N_GAPS_ACCEPTED.load(Relaxed)
    );
    eprintln!(
        "[predict] === split_by_ends routes: dense-cell={}  enclosure-fold={}  dropped={}  | folds WIDER than isoform extent={} ===",
        N_DENSE.load(Relaxed),
        N_ENCLOSED.load(Relaxed),
        N_DROP.load(Relaxed),
        N_ENCLOSED_WIDER.load(Relaxed)
    );
    let rejections = load(&N_REJ);
    let total = rejections.iter().sum::<u64>().max(1);
    eprintln!(
        "[predict] === place_on_backbone rejections: {} ===",
        (0..6)
            .map(|i| format!(
                "{}={} ({:.1}%)",
                REJECT_LABELS[i],
                rejections[i],
                100.0 * rejections[i] as f64 / total as f64
            ))
            .collect::<Vec<_>>()
            .join("  ")
    );
}
