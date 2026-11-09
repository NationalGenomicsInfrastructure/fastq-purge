//! Pipeline orchestration (M3: the exact pipeline; M4: the approx/bloom
//! pipeline; M5: the passes run on the global rayon pool sized by
//! `--threads` — see RUST_REWRITE.md §8, M5).

pub mod approx;
pub mod exact;

use crate::cli::{AppConfig, Method, ThreadingMethod};
use crate::logging::{Level, Logger};
use crate::membership::exact::ExactSet;

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    // One-time compatibility warning (RUST_REWRITE.md §7): the
    // --threading-method value is accepted and ignored; the Rust
    // implementation is thread-based and --threads sets the pool size.
    if config.threading_method != ThreadingMethod::MpPool {
        config.log.warning(&format!(
            "'--threading-method {}' accepted and ignored (the Rust \
             implementation is thread-based; '--threads {}' sets the pool \
             size)",
            config.threading_method.as_cli_name(),
            config.threads
        ));
    }

    match config.method {
        Method::Exact => exact::run(config),
        Method::Approx => approx::run(config),
    }
}

/// Log lines collected by a pass worker, in the order the sequential run
/// would have printed them (the M5 deterministic-log rule: workers never
/// print; the orchestrator flushes each worker's lines in canonical
/// lane/file order, so the log is byte-identical to the sequential run
/// for any `--threads`).
#[derive(Debug, Clone, Default)]
pub struct Lines {
    v: Vec<(Level, String)>,
}

impl Lines {
    pub fn debug(&mut self, msg: impl Into<String>) {
        self.v.push((Level::Debug, msg.into()));
    }

    pub fn info(&mut self, msg: impl Into<String>) {
        self.v.push((Level::Info, msg.into()));
    }

    pub fn warning(&mut self, msg: impl Into<String>) {
        self.v.push((Level::Warning, msg.into()));
    }

    pub fn into_vec(self) -> Vec<(Level, String)> {
        self.v
    }
}

/// One pass item's collected output: the log lines (Ok) or the lines
/// collected up to the failure plus the error (Err). The orchestrator
/// flushes these in canonical lane/file order (M5 deterministic-log
/// rule).
pub type Block = Result<Lines, (Lines, anyhow::Error)>;

/// Pass 2's per-file item: the per-lane partial dups sets plus the log
/// lines.
pub type ScanBlock = Result<(Vec<ExactSet>, Lines), (Lines, anyhow::Error)>;

/// Print collected worker lines through the level-gated logger.
pub(crate) fn flush_lines(log: &Logger, lines: &[(Level, String)]) {
    for (level, msg) in lines {
        log.emit(*level, msg);
    }
}
