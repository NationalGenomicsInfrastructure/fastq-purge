//! Pipeline orchestration (milestone M3: the exact pipeline, single
//! threaded; M4 adds the bloom/approx pipeline; M5 adds rayon).

pub mod exact;

use crate::cli::{AppConfig, Method, ThreadingMethod};

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    // One-time compatibility warning (RUST_REWRITE.md §7): the
    // --threading-method value is accepted and ignored, and M3's
    // pipeline is single threaded so --threads > 1 is not honoured yet.
    if config.threading_method != ThreadingMethod::MpPool || config.threads != 1 {
        config.log.warning(&format!(
            "'--threading-method {}' accepted and ignored (the Rust \
             implementation is thread-based); '--threads {}' ignored by \
             the single-threaded pipeline (multi-threading: M5)",
            config.threading_method.as_cli_name(),
            config.threads
        ));
    }

    match config.method {
        Method::Exact => exact::run(config),
        Method::Approx => Err(anyhow::anyhow!(
            "fastq-purge-rs M3: the approx (bloom) pipeline lands in M4; \
             the exact pipeline is complete"
        )),
    }
}
