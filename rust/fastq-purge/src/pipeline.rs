//! Pipeline orchestration.
//!
//! Milestone M1 (RUST_REWRITE.md §8) delivers CLI + discovery/validation
//! only; the passes of §5 land in M3 (exact) and M4 (approx).

use crate::cli::AppConfig;

pub fn run(_config: &AppConfig) -> anyhow::Result<()> {
    Err(anyhow::anyhow!(
        "fastq-purge-rs M1: pipeline not implemented yet (exact: M3, approx: M4); \
         discovery/validation complete"
    ))
}
