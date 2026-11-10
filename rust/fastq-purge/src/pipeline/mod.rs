//! Pipeline orchestration (M3: the exact pipeline; M4: the approx/bloom
//! pipeline; M5: the passes run on the global rayon pool sized by
//! `--threads` — see RUST_REWRITE.md §8, M5).

pub mod approx;
pub mod exact;
pub mod sheet;

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

    // M7: the SampleSheet input mode runs the branch's per-lane contract
    // instead of the canonical global 3-pass (D8 superset).
    if config.mode == crate::cli::RunMode::Sheet {
        return sheet::run(config);
    }

    match config.method {
        Method::Exact => exact::run(config),
        Method::Approx => approx::run(config),
    }
}

/// A pipeline error that has already been reported through the logger
/// (the branch-style ERROR line on a sheet-mode write failure). `main`
/// prints pipeline errors to stderr — except when the Display is empty,
/// which marks an already-logged error.
#[derive(Debug)]
pub struct Logged;

impl std::fmt::Display for Logged {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        // Intentionally empty: the message was already logged.
        f.write_str("")
    }
}

impl std::error::Error for Logged {}

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

    /// A copy of the buffer (for flushes that must not consume it).
    pub(crate) fn to_vec(&self) -> Vec<(Level, String)> {
        self.v.clone()
    }

    /// Append another block's lines (M7: combining a worker block with
    /// the lane's running block).
    pub(crate) fn append(&mut self, other: Lines) {
        self.v.extend(other.v);
    }

    /// Append the `memory_usage` line pair (M7 sheet mode; /proc/meminfo,
    /// live values — the differential normalizes these lines).
    pub(crate) fn append_mem(&mut self, pid: Option<u32>) {
        for msg in sheet::mem_lines(pid) {
            self.debug(msg);
        }
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
