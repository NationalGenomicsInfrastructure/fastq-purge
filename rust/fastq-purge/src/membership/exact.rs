//! Exact membership: a `HashSet<u64>` over xxhash64(name0).
//!
//! The 64-bit collision probability is quantified in RUST_REWRITE.md §6
//! (≈2.7e-4 chance of any colliding pair at 100 M names) and accepted for
//! "exact" mode (risk R2); an optional string-storing `--strict` mode is a
//! possible later addition, not built speculatively.

use std::collections::HashSet;

use xxhash_rust::xxh64::xxh64;

pub struct ExactSet {
    inner: HashSet<u64>,
}

impl ExactSet {
    pub fn new() -> Self {
        Self {
            inner: HashSet::new(),
        }
    }

    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            inner: HashSet::with_capacity(capacity),
        }
    }

    /// Insert the hash of a name0.
    pub fn insert_name(&mut self, name0: &[u8]) {
        self.inner.insert(xxh64(name0, 0));
    }

    /// True if the name0 was inserted.
    pub fn contains_name(&self, name0: &[u8]) -> bool {
        self.inner.contains(&xxh64(name0, 0))
    }

    pub fn len(&self) -> usize {
        self.inner.len()
    }

    pub fn is_empty(&self) -> bool {
        self.inner.is_empty()
    }
}

impl Default for ExactSet {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn insert_and_contains() {
        let mut s = ExactSet::new();
        assert!(s.is_empty());
        s.insert_name(b"@M00150:2:AABCCC:1:1107:1607:1082");
        assert!(s.contains_name(b"@M00150:2:AABCCC:1:1107:1607:1082"));
        assert!(!s.contains_name(b"@M00150:2:AABCCC:1:1107:1607:1083"));
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn duplicates_collapse() {
        let mut s = ExactSet::new();
        for _ in 0..1000 {
            s.insert_name(b"samename");
        }
        assert_eq!(s.len(), 1);
    }

    #[test]
    fn empty_name0_hashes_fine() {
        let mut s = ExactSet::new();
        s.insert_name(b"");
        assert!(s.contains_name(b""));
        assert!(!s.contains_name(b"x"));
    }

    #[test]
    fn distinct_name0s_collide_rarely() {
        // The set hashes whatever slice it is given; name0 *extraction*
        // (stop at first space) happens in the reader. Here: 10k distinct
        // name0s must all land in distinct buckets (64-bit space).
        let mut s = ExactSet::new();
        for i in 0..10_000u64 {
            s.insert_name(format!("M00150:2:AABCCC:1:1107:{i}:1082").as_bytes());
        }
        assert_eq!(s.len(), 10_000);
        for i in 0..10_000u64 {
            assert!(s.contains_name(format!("M00150:2:AABCCC:1:1107:{i}:1082").as_bytes()));
        }
        assert!(!s.contains_name(b"M00150:2:AABCCC:1:1107:99999:1082"));
    }
}
