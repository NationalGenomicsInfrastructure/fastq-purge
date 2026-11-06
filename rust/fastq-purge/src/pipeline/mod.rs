//! Pipeline orchestration (M3: the exact pipeline; M4: the approx/bloom
//! pipeline; both single threaded; M5 adds rayon).

pub mod approx;
pub mod exact;

use crate::cli::{AppConfig, Method, ThreadingMethod};

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    // One-time compatibility warning (RUST_REWRITE.md §7): the
    // --threading-method value is accepted and ignored, and the
    // pipelines are single threaded so --threads > 1 is not honoured yet.
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
        Method::Approx => approx::run(config),
    }
}
