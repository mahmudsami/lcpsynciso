// ─────────────────────────────────────────────────────────────────────────────
// DUPLICATED FROM synpact/src/hash.rs (verbatim, trimmed to what seedhist needs).
// Keep in sync with synpact so the seed/block hashes match exactly.
// ─────────────────────────────────────────────────────────────────────────────

/// The atom value for a window: base-4 encoding of the (forward) k-mer.
#[inline]
pub fn atom_value(atom: &[u8]) -> u64 {
    encode_smer(atom)
}

// Forward-only DNA rolling hash (rotation-XOR; NT-hash forward component).
pub const DNA_SEED: [u64; 256] = {
    let mut t = [0u64; 256];
    t[b'A' as usize] = 0x9e3779b97f4a7c15;
    t[b'a' as usize] = 0x9e3779b97f4a7c15;
    t[b'C' as usize] = 0x6c62272e07bb0142;
    t[b'c' as usize] = 0x6c62272e07bb0142;
    t[b'G' as usize] = 0xbf58476d1ce4e5b9;
    t[b'g' as usize] = 0xbf58476d1ce4e5b9;
    t[b'T' as usize] = 0x94d049bb133111eb;
    t[b't' as usize] = 0x94d049bb133111eb;
    t
};

/// Rolling forward-only DNA hash over a window of `w` consecutive bases.
pub struct DnaHashFwd<'a> {
    pub seq:  &'a [u8],
    pub w:    usize,
    pub h:    u64,
    pub n_in_window: u32,
    pub pos:  usize,
}

impl<'a> DnaHashFwd<'a> {
    pub fn new(seq: &'a [u8], w: usize) -> Option<Self> {
        if seq.len() < w { return None; }
        let mut h = 0u64;
        let mut n_count = 0u32;
        for &b in &seq[..w] {
            if DNA_SEED[b as usize] == 0 { n_count += 1; }
            h = h.rotate_left(1) ^ DNA_SEED[b as usize];
        }
        Some(DnaHashFwd { seq, w, h, n_in_window: n_count, pos: 0 })
    }
}

impl<'a> Iterator for DnaHashFwd<'a> {
    type Item = u64;
    #[inline]
    fn next(&mut self) -> Option<u64> {
        if self.pos + self.w > self.seq.len() { return None; }
        let h = if self.n_in_window > 0 { u64::MAX } else { self.h };
        if self.pos + self.w < self.seq.len() {
            let out = self.seq[self.pos];
            let in_ = self.seq[self.pos + self.w];
            if DNA_SEED[out as usize] == 0 { self.n_in_window -= 1; }
            if DNA_SEED[in_  as usize] == 0 { self.n_in_window += 1; }
            self.h = self.h.rotate_left(1)
                ^ DNA_SEED[out as usize].rotate_left(self.w as u32)
                ^ DNA_SEED[in_  as usize];
        }
        self.pos += 1;
        Some(h)
    }
}

#[inline]
pub fn base4(b: u8) -> u64 {
    match b.to_ascii_uppercase() {
        b'A' => 0, b'C' => 1, b'G' => 2, b'T' => 3, _ => 0,
    }
}

#[inline]
pub fn encode_smer(smer: &[u8]) -> u64 {
    smer.iter().fold(0u64, |acc, &b| acc * 4 + base4(b))
}

pub const LEVEL_DOMAINS: [u64; 8] = [
    0x9e3779b97f4a7c15, // L1
    0x6c62272e07bb0142, // L2
    0xd2a98b26625eee7b, // L3
    0xa3b195354a2b7623, // L4
    0x1b03738712fad5c9, // L5
    0xc4ceb9fe1a85ec53, // L6
    0x517cc1b727220a95, // L7
    0x9be1aa58ba6a9f81, // L8
];

#[inline]
pub fn level_domain(level_0idx: usize) -> u64 {
    LEVEL_DOMAINS.get(level_0idx).copied().unwrap_or_else(||
        LEVEL_DOMAINS[7].wrapping_add(
            (level_0idx as u64).wrapping_mul(0x6c62272e07bb0142)))
}

#[inline]
pub fn block_hash_iter_domain(len: usize, vals: impl Iterator<Item = u64>, domain: u64) -> u64 {
    let mut h: u64 = (len as u64)
        .wrapping_mul(0xbf58476d1ce4e5b9)
        ^ domain;
    for v in vals {
        h = h.rotate_left(13).wrapping_add(v);
        h ^= h >> 30;
        h = h.wrapping_mul(0xbf58476d1ce4e5b9);
        h ^= h >> 27;
        h = h.wrapping_mul(0x94d049bb133111eb);
        h ^= h >> 31;
    }
    h
}

#[inline]
pub fn block_hash_for_level(vals: &[u64], level_0idx: usize) -> u64 {
    block_hash_iter_domain(vals.len(), vals.iter().copied(), level_domain(level_0idx))
}
