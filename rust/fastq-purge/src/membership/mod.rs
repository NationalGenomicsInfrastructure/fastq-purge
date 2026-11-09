//! Membership structures for the pipeline passes.
//!
//! M2: the exact set (`exact`). M4: the bloom filter (`bloom`).

pub mod bloom;
pub mod exact;

pub use bloom::Bloom64;
pub use exact::ExactSet;
