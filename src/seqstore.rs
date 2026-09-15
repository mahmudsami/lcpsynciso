//! Compact in-RAM read store: 2-bit-packed DNA plus read names.
//!
//! Lets the fused `predict` pipeline read the input file exactly once, hold every
//! read in memory, and hand sequences to isoform resolution without ever
//! re-reading or re-extracting them. ACGT pack to 2 bits/base (~4x smaller than
//! raw bytes); the rare non-ACGT byte (e.g. `N`) is packed as `A` and its
//! original value recorded in a sparse exception list, so `get` restores the
//! read verbatim.

use std::collections::HashMap;

/// 2-bit code for a base, or `None` for non-ACGT.
#[inline]
fn code(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn base(c: u64) -> u8 {
    match c & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

#[derive(Default)]
pub struct SeqStore {
    names: Vec<String>,
    packed: Vec<u64>,   // 32 bases per word, all reads concatenated
    starts: Vec<u64>,   // base offset of read i into the global packing
    lengths: Vec<u32>,  // read length in bases
    total_bases: u64,   // running base count = next write position
    exceptions: HashMap<u32, Vec<(u32, u8)>>, // read idx -> (pos, original byte) for non-ACGT
}

impl SeqStore {
    pub fn with_capacity(n_reads: usize, total_bases: u64) -> Self {
        SeqStore {
            names: Vec::with_capacity(n_reads),
            packed: Vec::with_capacity((total_bases as usize / 32) + 1),
            starts: Vec::with_capacity(n_reads),
            lengths: Vec::with_capacity(n_reads),
            total_bases: 0,
            exceptions: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn name(&self, i: usize) -> &str {
        &self.names[i]
    }

    pub fn read_len(&self, i: usize) -> usize {
        self.lengths[i] as usize
    }

    /// Append a read; returns its index.
    pub fn push(&mut self, name: String, seq: &[u8]) -> u32 {
        let idx = self.names.len() as u32;
        let start = self.total_bases;
        self.starts.push(start);
        self.lengths.push(seq.len() as u32);
        self.names.push(name);

        let mut exc: Vec<(u32, u8)> = Vec::new();
        for (p, &b) in seq.iter().enumerate() {
            let global = self.total_bases + p as u64;
            let word = (global / 32) as usize;
            let shift = ((global % 32) * 2) as u64;
            if word >= self.packed.len() {
                self.packed.push(0);
            }
            let c = match code(b) {
                Some(c) => c,
                None => {
                    exc.push((p as u32, b));
                    0
                }
            };
            self.packed[word] |= c << shift;
        }
        self.total_bases += seq.len() as u64;
        if !exc.is_empty() {
            self.exceptions.insert(idx, exc);
        }
        idx
    }

    /// Reconstruct read `i`'s sequence bytes (upper-case ACGT plus restored
    /// non-ACGT bytes).
    pub fn get(&self, i: usize) -> Vec<u8> {
        let start = self.starts[i];
        let len = self.lengths[i] as usize;
        let mut out = Vec::with_capacity(len);
        for p in 0..len {
            let global = start + p as u64;
            let word = self.packed[(global / 32) as usize];
            let shift = ((global % 32) * 2) as u64;
            out.push(base(word >> shift));
        }
        if let Some(exc) = self.exceptions.get(&(i as u32)) {
            for &(p, b) in exc {
                out[p as usize] = b;
            }
        }
        out
    }
}
