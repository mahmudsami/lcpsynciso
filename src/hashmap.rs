//! Identity-hashed `HashMap<u64, V>` used on the hot seed-counting paths.
//! Seed hashes are already well-mixed, so re-hashing them is wasted work.

use std::collections::HashMap;
use std::hash::{BuildHasherDefault, Hasher};

#[derive(Default)]
pub struct IdHasher(u64);

impl Hasher for IdHasher {
    #[inline]
    fn finish(&self) -> u64 {
        self.0
    }
    #[inline]
    fn write(&mut self, _: &[u8]) {
        unreachable!("only u64 keys")
    }
    #[inline]
    fn write_u64(&mut self, n: u64) {
        self.0 = n;
    }
}

pub type FastMap<V> = HashMap<u64, V, BuildHasherDefault<IdHasher>>;

pub fn fmap<V>() -> FastMap<V> {
    FastMap::default()
}
