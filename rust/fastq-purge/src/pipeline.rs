//! Pipeline orchestration.
//!
//! Milestone M2 (RUST_REWRITE.md §8) delivers the FASTQ streaming reader/
//! writer and the exact membership set on top of M1's
//! CLI + discovery/validation; the passes of §5 land in M3 (exact) and
//! M4 (approx).

use crate::cli::AppConfig;

pub fn run(_config: &AppConfig) -> anyhow::Result<()> {
    Err(anyhow::anyhow!(
        "fastq-purge-rs M2: FASTQ streaming + name0 hashing complete \
         (fastq/, membership/); pipeline not implemented yet (exact: M3, \
         approx: M4)"
    ))
}
