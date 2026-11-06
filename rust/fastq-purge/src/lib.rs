//! fastq-purge (Rust rewrite). See `RUST_REWRITE.md` in the repository root
//! for the design, the Python-behaviour porting contract, and milestone
//! status.
//!
//! Crate layout (M2): `cli` + `discovery` (M1), `fastq` streaming reader/
//! writer + name0 handling (M2), `membership` (M2: exact set; M4: bloom),
//! `pipeline` (M3/M4). The `fastq-dump` bin is a dev/verification tool used
//! for differential testing against the Python implementation.

pub mod cli;
pub mod discovery;
pub mod fastq;
pub mod logging;
pub mod membership;
pub mod pipeline;
