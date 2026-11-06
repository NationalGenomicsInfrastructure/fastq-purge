//! Membership structures for the pipeline passes.
//!
//! M2: the exact set (`exact`). M4 adds the bloom filter (`bloom`).

pub mod exact;

pub use exact::ExactSet;
