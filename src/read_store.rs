//! Compact in-RAM read store: 2-bit-packed DNA plus read names.
//!
//! Lets the fused `predict` pipeline read the input file exactly once, hold every
//! read in memory, and hand sequences to isoform detection without ever
//! re-reading or re-extracting them. ACGT pack to 2 bits/base (~4x smaller than
//! raw bytes); the rare non-ACGT byte (e.g. `N`) is packed as `A` and its
//! original value recorded in a sparse non-ACGT list, so `read_seq` restores the
//! read verbatim.

use std::collections::HashMap;

/// 2-bit code for a base, or `None` for non-ACGT.
#[inline]
fn encode_base(b: u8) -> Option<u64> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

#[inline]
fn decode_base(c: u64) -> u8 {
    match c & 3 {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        _ => b'T',
    }
}

#[derive(Default)]
pub struct ReadStore {
    names: Vec<String>,
    packed: Vec<u64>,       // 32 bases per word, all reads concatenated
    base_offsets: Vec<u64>, // base offset of each read into the global packing
    lengths: Vec<u32>,      // read length in bases
    total_bases: u64,       // running base count = next write position
    non_acgt: HashMap<u32, Vec<(u32, u8)>>, // read idx -> (pos, original byte) for non-ACGT
}

impl ReadStore {
    pub fn with_capacity(expected_reads: usize, expected_bases: u64) -> Self {
        ReadStore {
            names: Vec::with_capacity(expected_reads),
            packed: Vec::with_capacity((expected_bases as usize / 32) + 1),
            base_offsets: Vec::with_capacity(expected_reads),
            lengths: Vec::with_capacity(expected_reads),
            total_bases: 0,
            non_acgt: HashMap::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.names.len()
    }

    pub fn is_empty(&self) -> bool {
        self.names.is_empty()
    }

    pub fn read_name(&self, read_idx: usize) -> &str {
        &self.names[read_idx]
    }

    pub fn read_len(&self, read_idx: usize) -> usize {
        self.lengths[read_idx] as usize
    }

    /// Append a read; returns its index.
    pub fn push(&mut self, name: String, seq: &[u8]) -> u32 {
        let idx = self.names.len() as u32;
        let start = self.total_bases;
        self.base_offsets.push(start);
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
            let c = match encode_base(b) {
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
            self.non_acgt.insert(idx, exc);
        }
        idx
    }

    /// Reconstruct read `read_idx`'s sequence bytes (upper-case ACGT plus restored
    /// non-ACGT bytes).
    pub fn read_seq(&self, read_idx: usize) -> Vec<u8> {
        let start = self.base_offsets[read_idx];
        let len = self.lengths[read_idx] as usize;
        let mut out = Vec::with_capacity(len);
        for p in 0..len {
            let global = start + p as u64;
            let word = self.packed[(global / 32) as usize];
            let shift = ((global % 32) * 2) as u64;
            out.push(decode_base(word >> shift));
        }
        if let Some(exc) = self.non_acgt.get(&(read_idx as u32)) {
            for &(p, b) in exc {
                out[p as usize] = b;
            }
        }
        out
    }
}
