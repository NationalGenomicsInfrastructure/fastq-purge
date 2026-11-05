//! Minimal level-gated logger mirroring the Python CLI's output contract.
//!
//! The Python tool logs through `RichHandler` (see `main.py`'s
//! `logging.basicConfig`), which renders every record on **stdout** as
//! `HH:MM:SS LEVEL message  main.py:LINE` with rich's wrapping and colour
//! (the `format="%(message)s"` argument is ignored by RichHandler).
//!
//! We reproduce the *message* contract — same messages, same order, same
//! `--log-level` gating, on stdout — but as plain bare lines. The rich
//! rendering (wall-clock timestamps, source locations, wrapping, colour)
//! is deliberately not reproduced: it is display chrome, it is
//! non-deterministic (timestamps), and it references Python source lines
//! that do not exist in the Rust rewrite.
//!
//! Fatal errors that abort the CLI (invalid path, M1 pipeline bail) also
//! use this logger, matching Python where those messages are logged
//! before `exit(1)`; only errors with no Python equivalent go to stderr.

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Level {
    Debug = 0,
    Info = 1,
    Warning = 2,
    Error = 3,
    Critical = 4,
}

impl Level {
    /// Parse the `--log-level` value. The Python CLI accepts exactly the
    /// five uppercase names (argparse choices are case-sensitive).
    pub fn parse(s: &str) -> Option<Level> {
        Some(match s {
            "DEBUG" => Level::Debug,
            "INFO" => Level::Info,
            "WARNING" => Level::Warning,
            "ERROR" => Level::Error,
            "CRITICAL" => Level::Critical,
            _ => return None,
        })
    }
}

#[derive(Debug, Clone)]
pub struct Logger {
    min: Level,
}

impl Logger {
    pub fn new(min: Level) -> Self {
        Self { min }
    }

    fn log(&self, level: Level, msg: &str) {
        if level >= self.min {
            println!("{msg}");
        }
    }

    pub fn debug(&self, msg: &str) {
        self.log(Level::Debug, msg);
    }

    pub fn info(&self, msg: &str) {
        self.log(Level::Info, msg);
    }

    pub fn warning(&self, msg: &str) {
        self.log(Level::Warning, msg);
    }

    pub fn error(&self, msg: &str) {
        self.log(Level::Error, msg);
    }

    pub fn critical(&self, msg: &str) {
        self.log(Level::Critical, msg);
    }
}
