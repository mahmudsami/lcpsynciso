//! Sequence comparison with block-aligner (SIMD, adaptive block size).
//!
//!   * [`flank`]: how much of a read's terminal flank the reference accounts for. An x-drop
//!     extension stops where the sequences stop agreeing, so the unaligned rest of the read
//!     is the unmatched flank, and whether the reference ran out tells an overhang from a
//!     divergence.
//!   * [`gap_runs`]: global alignment of an interior gap pinned between two anchors, read as
//!     substitution identity and longest indel.
//!
//! Buffers are thread-local and reused, since these run millions of times.

use std::cell::RefCell;

use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;

/// Match +1, mismatch -1. The match score must be positive, or an x-drop extension could
/// never gain score and would stop immediately.
static SW: NucMatrix = NucMatrix::new_simple(1, -1);
/// Affine gaps; block-aligner requires opening a gap to cost more than extending it.
const GAPS: Gaps = Gaps { open: -2, extend: -1 };

/// Block size range for the adaptive aligner.
const BS_LO: usize = 32;
const BS_HI: usize = 64;

/// Padded query and reference buffers, grown on demand.
struct Bufs {
    q: PaddedBytes,
    r: PaddedBytes,
    cap_q: usize,
    cap_r: usize,
}

impl Bufs {
    fn new() -> Self {
        Bufs {
            q: PaddedBytes::new::<NucMatrix>(BS_HI, BS_HI),
            r: PaddedBytes::new::<NucMatrix>(BS_HI, BS_HI),
            cap_q: BS_HI,
            cap_r: BS_HI,
        }
    }

    /// Make room for a query of `lq` and a reference of `lr` bases.
    fn fit(&mut self, lq: usize, lr: usize) {
        if lq > self.cap_q {
            self.cap_q = (lq * 2).max(BS_HI);
            self.q = PaddedBytes::new::<NucMatrix>(self.cap_q, BS_HI);
        }
        if lr > self.cap_r {
            self.cap_r = (lr * 2).max(BS_HI);
            self.r = PaddedBytes::new::<NucMatrix>(self.cap_r, BS_HI);
        }
    }
}

thread_local! {
    static BUFS: RefCell<Bufs> = RefCell::new(Bufs::new());
    /// Traceback buffer and its capacity, grown on demand.
    static CIGAR: RefCell<(Cigar, usize)> = RefCell::new((Cigar::new(BS_HI, BS_HI), BS_HI));
}

/// Result of extending a read flank into the reference flank.
pub struct FlankFit {
    /// Read-flank bases the reference does not account for.
    pub unmatched: usize,
    /// Reference bases consumed by the alignment.
    pub consumed: usize,
    /// The reference flank ran out: the read extends past it (an overhang). False means the
    /// alignment stopped while reference sequence remained (the sequences diverge).
    pub reference_exhausted: bool,
}

/// Extend read flank `read` into reference flank `reference`, outward from their anchor. With
/// `rev`, both are aligned reversed, so a 5' flank also starts at the anchor.
pub fn flank(read: &[u8], reference: &[u8], rev: bool) -> FlankFit {
    if read.is_empty() {
        return FlankFit { unmatched: 0, consumed: 0, reference_exhausted: reference.is_empty() };
    }
    if reference.is_empty() {
        return FlankFit { unmatched: read.len(), consumed: 0, reference_exhausted: true };
    }
    BUFS.with(|b| {
        let b = &mut *b.borrow_mut();
        b.fit(read.len(), reference.len());
        if rev {
            b.q.set_bytes_rev::<NucMatrix>(read, BS_HI);
            b.r.set_bytes_rev::<NucMatrix>(reference, BS_HI);
        } else {
            b.q.set_bytes::<NucMatrix>(read, BS_HI);
            b.r.set_bytes::<NucMatrix>(reference, BS_HI);
        }
        // X-drop scaled to the flank: too small cuts a real match short, too large walks
        // through divergent sequence.
        let xdrop = ((read.len() / 4) as i32).max(20);
        let mut blk = Block::<false, true>::new(read.len(), reference.len(), BS_HI);
        blk.align(&b.q, &b.r, &SW, GAPS, BS_LO..=BS_HI, xdrop);
        let res = blk.res();
        FlankFit {
            unmatched: read.len().saturating_sub(res.query_idx),
            consumed: res.reference_idx.min(reference.len()),
            reference_exhausted: res.reference_idx >= reference.len(),
        }
    })
}

/// Global alignment of `q` against `r` with traceback. `read_trace` gets the trace, the
/// padded query and reference, the alignment result, and the reusable CIGAR buffer.
fn traced_global<T>(
    q: &[u8],
    r: &[u8],
    read_trace: impl FnOnce(&Trace, &PaddedBytes, &PaddedBytes, AlignResult, &mut Cigar) -> T,
) -> T {
    BUFS.with(|bf| {
        CIGAR.with(|cg| {
            let bf = &mut *bf.borrow_mut();
            let (cigar, cap) = &mut *cg.borrow_mut();
            bf.fit(q.len(), r.len());
            bf.q.set_bytes::<NucMatrix>(q, BS_HI);
            bf.r.set_bytes::<NucMatrix>(r, BS_HI);
            // `Cigar::new` sizes for query + reference, so grow it when needed.
            let need = q.len() + r.len();
            if need > *cap {
                *cap = need * 2;
                *cigar = Cigar::new(*cap, *cap);
            }
            let mut blk = Block::<true, false>::new(q.len(), r.len(), BS_HI);
            blk.align(&bf.q, &bf.r, &SW, GAPS, BS_LO..=BS_HI, 0);
            let res = blk.res();
            read_trace(blk.trace(), &bf.q, &bf.r, res, cigar)
        })
    })
}

/// Summary of a global alignment between two interior gaps.
pub struct GapStats {
    /// Longest single run of inserted or deleted bases.
    pub max_indel: usize,
    /// Matching columns.
    pub eq: usize,
    /// Mismatching columns.
    pub mismatch: usize,
}

impl GapStats {
    /// Identity over aligned columns only, Eq / (Eq + X). Indel columns are excluded, so a
    /// length difference does not lower it; `max_indel` reports that separately.
    pub fn ident(&self) -> f64 {
        let d = self.eq + self.mismatch;
        if d == 0 { 1.0 } else { self.eq as f64 / d as f64 }
    }
}

/// Align two interior gaps globally and summarise the alignment.
pub fn gap_runs(a: &[u8], b: &[u8]) -> GapStats {
    if a.is_empty() || b.is_empty() {
        return GapStats { max_indel: a.len().max(b.len()), eq: 0, mismatch: 0 };
    }
    traced_global(a, b, |trace, q, r, res, cigar| {
        // `cigar_eq` separates matches (=) from mismatches (X); plain `cigar` merges them as M.
        trace.cigar_eq(q, r, res.query_idx, res.reference_idx, cigar);
        let (mut max_indel, mut eq, mut mismatch) = (0, 0, 0);
        for i in 0..cigar.len() {
            let op = cigar.get(i);
            match op.op {
                Operation::Eq => eq += op.len,
                Operation::X => mismatch += op.len,
                Operation::I | Operation::D => max_indel = max_indel.max(op.len),
                _ => {}
            }
        }
        GapStats { max_indel, eq, mismatch }
    })
}
