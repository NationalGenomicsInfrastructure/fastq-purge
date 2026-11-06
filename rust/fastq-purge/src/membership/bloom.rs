//! Approximate membership: a bloom filter over name0s (milestone M4).
//!
//! Sizing (RUST_REWRITE.md §6, corrected in §8/M4 decision D4):
//!
//! ```text
//! m = ⌈1.2 · n · (−ln fpr) / (ln 2)²⌉ bits     (rbloom-parity formula,
//!                                                 +20% margin)
//! k = max(1, round(m / n · ln 2))                (optimal k for the size)
//! ```
//!
//! The design draft fixed k = 2 (the two halves of xxh128); at 19.17
//! bits/item that combination yields a ≈1% structural FPR, not the
//! promised `--fpr`. We keep the *size* formula (so the memory table in
//! §6 holds: 19.17 bits/item at fpr 1e-4) and use the optimal k with
//! Kirsch–Mitzenmacher double hashing from two 64-bit xxhashes:
//!
//! ```text
//! h1 = xxh64(name0, 0),  h2 = xxh64(name0, 0x9E3779B97F4A7C15)
//! probe_i = (h1 + i · h2) mod m,  i = 0..k
//! ```
//!
//! Achieved structural FPR ≈ (1 − e^(−kn/m))^k, which with the 20% margin
//! is *better* than `--fpr`; the empirical check is the §8 property test.
//! The 2024 double-hash quality caveat (Ebeling et al.) is not material
//! here: probes only gate *detection* of duplicates in pass 2 of the
//! canonical pipeline (decision D4), where the duplicate set itself is
//! exact, so structural false positives cannot change which reads are
//! removed.

use xxhash_rust::xxh64::xxh64;

/// Second xxh64 seed for the double hash (golden-ratio constant).
const SEED2: u64 = 0x9E37_79B9_7F4A_7C15;

/// Both 64-bit halves of the double hash for a name0, computed once and
/// shared across several filters (pass 2 tests one name0 against all K
/// lane filters).
pub fn hash_name0(name0: &[u8]) -> (u64, u64) {
    (xxh64(name0, 0), xxh64(name0, SEED2))
}

pub struct Bloom64 {
    /// Bit array, word-major (bit `b` lives in `words[b / 64]` bit `b % 64`).
    words: Vec<u64>,
    /// Total number of bits (`words.len() * 64`; the tail of the last word
    /// is unused when `m` is not a multiple of 64).
    m: u64,
    /// Number of probes.
    k: u32,
    /// Number of insert calls (diagnostic; duplicates included).
    observed: u64,
    /// `n_items` from the constructor (diagnostic).
    capacity: u64,
}

impl Bloom64 {
    /// Size `m` and `k` per the module doc for `n_items` expected items and
    /// target FPR `fpr`. `n_items` is floored at 1; `fpr` is clamped to
    /// `[1e-9, 0.99]` (argparse, like rbloom, does not validate the flag).
    pub fn new(n_items: u64, fpr: f64) -> Self {
        let n = n_items.max(1);
        let fpr = fpr.clamp(1e-9, 0.99);
        let ln2 = std::f64::consts::LN_2;
        let m = ((n as f64) * 1.2 * (-fpr.ln()) / (ln2 * ln2)).ceil() as u64;
        let k = ((m as f64 / n as f64) * ln2).round() as u32;
        let k = k.max(1);
        let words = vec![0u64; (m.div_ceil(64)) as usize];
        Self {
            words,
            m,
            k,
            observed: 0,
            capacity: n,
        }
    }

    fn set_bit(&mut self, b: u64) {
        let w = b / 64;
        let bit = 1u64 << (b % 64);
        self.words[w as usize] |= bit;
    }

    fn get_bit(&self, b: u64) -> bool {
        let w = b / 64;
        let bit = 1u64 << (b % 64);
        self.words[w as usize] & bit != 0
    }

    /// Insert the precomputed double hash.
    pub fn insert_h(&mut self, h1: u64, h2: u64) {
        self.observed += 1;
        for i in 0..self.k {
            let probe = h1.wrapping_add(h2.wrapping_mul(i as u64));
            self.set_bit(probe % self.m);
        }
    }

    /// Membership test on a precomputed double hash (avoids re-hashing when
    /// one name0 is tested against several filters).
    pub fn contains_h(&self, h1: u64, h2: u64) -> bool {
        for i in 0..self.k {
            let probe = h1.wrapping_add(h2.wrapping_mul(i as u64));
            if !self.get_bit(probe % self.m) {
                return false;
            }
        }
        true
    }

    /// Insert the hash of a name0.
    pub fn insert_name(&mut self, name0: &[u8]) {
        let (h1, h2) = hash_name0(name0);
        self.insert_h(h1, h2);
    }

    /// True if the name0 was inserted.
    pub fn contains_name(&self, name0: &[u8]) -> bool {
        let (h1, h2) = hash_name0(name0);
        self.contains_h(h1, h2)
    }

    /// Total bits allocated.
    pub fn bits(&self) -> u64 {
        self.m
    }

    /// Number of probes (optimal k for the allocated size).
    pub fn k(&self) -> u32 {
        self.k
    }

    /// `n_items` as passed to `new` (diagnostic; see RUST_REWRITE.md R5).
    pub fn capacity(&self) -> u64 {
        self.capacity
    }

    /// Number of insert calls so far (diagnostic).
    pub fn observed(&self) -> u64 {
        self.observed
    }

    /// Allocated bytes (whole words).
    pub fn estimated_bytes(&self) -> usize {
        self.words.len() * 8
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sizing_matches_formulas() {
        // n = 1e6, fpr = 1e-4: m = ceil(1.2e6 × 19.1699) = 23,004,141 bits,
        // k = round(23.0041 × ln2) = 16.
        let b = Bloom64::new(1_000_000, 1e-4);
        assert_eq!(b.bits(), 23_004_141);
        assert_eq!(b.k(), 16);
        assert_eq!(b.estimated_bytes(), 2_875_520);
    }

    #[test]
    fn sizing_degenerates_safely() {
        let b = Bloom64::new(0, 1e-4); // n floored at 1
        assert!(b.bits() > 0 && b.k() >= 1);
        let b = Bloom64::new(10, 5.0); // fpr clamped to 0.99
        assert!(b.bits() > 0 && b.k() >= 1);
    }

    #[test]
    fn no_false_negatives() {
        let mut b = Bloom64::new(10_000, 1e-4);
        for i in 0..10_000u64 {
            b.insert_name(format!("name-{i}").as_bytes());
        }
        for i in 0..10_000u64 {
            assert!(b.contains_name(format!("name-{i}").as_bytes()));
        }
    }

    #[test]
    fn empirical_fpr_within_target() {
        // §8 property test (data-structure level): insert 100k, probe 100k
        // non-members, measured FPR ≤ target (the 20% margin makes the
        // achieved FPR ~5× better than the target).
        let n = 100_000u64;
        let target = 1e-4;
        let mut b = Bloom64::new(n, target);
        for i in 0..n {
            b.insert_name(format!("member-{i}").as_bytes());
        }
        let mut fp = 0u64;
        for i in 0..n {
            if b.contains_name(format!("nonmember-{i}").as_bytes()) {
                fp += 1;
            }
        }
        let measured = fp as f64 / n as f64;
        assert!(
            measured <= target,
            "measured FPR {measured:.6} exceeds target {target}"
        );
    }
}
