//! CLI definitions and validation entry point (milestone M1).
//!
//! Flags, defaults and choices mirror `fastq_purge.main.parse_args()`
//! (Python, argparse) exactly; discovery/validation behaviour mirrors
//! `validate_args()`. See RUST_REWRITE.md §7 for the parity checklist.

use std::path::{Path, PathBuf};

use clap::Parser;

use crate::discovery::{self, UndeterminedLane};
use crate::logging::{Level, Logger};

/// The `--method` choice (parity with argparse `choices`).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "lower")]
pub enum Method {
    Exact,
    Approx,
}

/// The `--threading-method` choice.
///
/// The Python flag selects a process-pool backend; the Rust implementation
/// uses threads only, so the value is accepted and ignored (one-time
/// warning, per RUST_REWRITE.md §7).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "kebab-case")]
pub enum ThreadingMethod {
    Loky,
    MpPool,
    MpManager,
}

impl ThreadingMethod {
    /// The name as it appears on the command line.
    ///
    /// Used by the pipeline phase (M3/M4) for the one-time "value
    /// accepted and ignored" warning; unused in M1.
    #[allow(dead_code)]
    pub fn as_cli_name(self) -> &'static str {
        match self {
            ThreadingMethod::Loky => "loky",
            ThreadingMethod::MpPool => "mp_pool",
            ThreadingMethod::MpManager => "mp_manager",
        }
    }
}

/// `--log-level` choice (parsed case-sensitively, like argparse choices).
#[derive(Debug, Clone, Copy, PartialEq, Eq, clap::ValueEnum)]
#[clap(rename_all = "UPPER")]
pub enum LogLevel {
    Debug,
    Info,
    Warning,
    Error,
    Critical,
}

impl LogLevel {
    pub fn level(self) -> Level {
        match self {
            LogLevel::Debug => Level::Debug,
            LogLevel::Info => Level::Info,
            LogLevel::Warning => Level::Warning,
            LogLevel::Error => Level::Error,
            LogLevel::Critical => Level::Critical,
        }
    }
}

#[derive(Parser, Debug)]
#[command(
    name = "fastq-purge-rs",
    version,
    about = "Package to purge fastq files from unwanted reads",
    override_usage = "fastq-purge-rs [options]"
)]
pub struct Cli {
    /// Path to the target fastq file(s) to purge. It can be gzipped or not.
    /// It can be a single file or a directory. If a directory is provided,
    /// all files in the directory that match the pattern '*.fq*' or
    /// '*.fastq*' will be used as target files. The search is non-recursive.
    #[arg(long, required = true)]
    pub undetermined_path: PathBuf,

    /// Path to the output fastq file
    #[arg(long)]
    pub output_path: Option<String>,

    /// Path to the bloom filter sources file(s). They can be gzipped or
    /// not. It can be a single file or multiple files separated by spaces,
    /// or a directory. If a directory is provided, all files in the
    /// directory that match the pattern '*.fq*' or '*.fastq*' will be used
    /// as bloom filter sources. The search is recursive.
    #[arg(long, required = true, num_args = 1..)]
    pub assigned_path: Vec<PathBuf>,

    /// Method to use for filtering. The 'exact' method uses a set to store
    /// the undetermined reads. The 'approx' method uses a bloom filter.
    /// The 'exact' method is more resource intensive, but it is more
    /// accurate. The 'approx' method is less resource intensive, but it has
    /// a false positive rate associated with it.
    #[arg(long, value_enum, default_value_t = Method::Exact)]
    pub method: Method,

    /// The estimated number of items in the bloom filter. It is used
    /// together with the false positive rate to calculate the size of the
    /// bloom filter. The default is 1 billion. This is ignored if the
    /// method is 'exact'.
    #[arg(long, default_value_t = 1_000_000_000)]
    pub max_items: u64,

    /// False positive rate for the bloom filter. The default is 0.0001.
    #[arg(long, default_value_t = 0.0001)]
    pub fpr: f64,

    /// Threading library and method to use for processing (used for testing
    /// purposes). Accepted for script compatibility; the Rust
    /// implementation is thread-based and ignores this value.
    #[arg(long, value_enum, default_value_t = ThreadingMethod::MpPool)]
    pub threading_method: ThreadingMethod,

    /// Number of threads to use for processing
    #[arg(long, default_value_t = 1)]
    pub threads: u64,

    /// Set the logging level
    #[arg(long, value_enum, default_value_t = LogLevel::Info)]
    pub log_level: LogLevel,
}

/// Validated, ready-to-run configuration (mirrors the Python
/// `args` namespace after `validate_args`).
///
/// The pipeline fields are consumed from M3/M4 onward; M1's stub
/// pipeline reads none of them.
#[allow(dead_code)]
pub struct AppConfig {
    pub log: Logger,
    /// lane key (e.g. "L001") -> (R1 or the single file, R2 or None), in
    /// discovery order (parity with the Python dict insertion order).
    pub undetermined: Vec<(String, UndeterminedLane)>,
    /// lane key -> (output R1, output R2 or None), same order.
    pub output: Vec<(String, (PathBuf, Option<PathBuf>))>,
    /// Directory the output files live in (the resolved `--output-path`).
    pub output_dir: PathBuf,
    /// lane key -> assigned source files, in discovery order.
    pub assigned: Vec<(String, Vec<PathBuf>)>,
    pub method: Method,
    pub max_items: u64,
    pub fpr: f64,
    pub threads: u64,
    pub threading_method: ThreadingMethod,
}

/// Validate `Cli` into an `AppConfig`, emitting the same log lines (same
/// content, same order) as the Python `validate_args`.
///
/// Returns `Err(msg, exit_code, already_logged)`: errors with
/// `already_logged == true` were already emitted to stdout by the logger
/// (parity with the Python CLI, which logs the "neither file nor
/// directory" message before exiting). argparse-style usage errors exit 2
/// before this point, handled by clap.
/// `(message, exit_code 1, not yet logged)`.
fn err(msg: String) -> (String, i32, bool) {
    (msg, 1, false)
}

/// `(message, exit_code 1, already logged to stdout)`.
fn logged_err(msg: String) -> (String, i32, bool) {
    (msg, 1, true)
}

pub fn validate(cli: Cli) -> Result<AppConfig, (String, i32, bool)> {
    let log = Logger::new(cli.log_level.level());

    // ---- undetermined -------------------------------------------------
    let un_files =
        discovery::explode_path(&cli.undetermined_path, false, &log).map_err(logged_err)?;

    let undetermined = discovery::create_paired_dict(un_files).map_err(err)?;
    if undetermined.is_empty() {
        // The Python CLI raises IndexError here (empty dict -> keys()[0]);
        // we fail cleanly instead.
        return Err(err(format!(
            "No undetermined files found in '{}' (pattern '*.fq*' or '*.fastq*')",
            cli.undetermined_path.display()
        )));
    }

    log.info("Undetermined files:");
    for (lane, (f1, f2)) in &undetermined {
        log.info(&format!("    {lane}: '{}'", f1.display()));
        match f2 {
            Some(f) => log.info(&format!("    {lane}: '{}'", f.display())),
            None => log.info(&format!("    {lane}: 'None'")),
        }
    }

    // ---- output path ---------------------------------------------------
    let output_dir = match &cli.output_path {
        Some(p) => PathBuf::from(p),
        None => undetermined
            .first()
            .expect("checked non-empty above")
            .1
             .0
            .parent()
            .unwrap_or(Path::new("."))
            .to_path_buf(),
    };
    if !output_dir.is_dir() {
        log.debug("The output path does not exist! It will be created...");
        std::fs::create_dir_all(&output_dir).map_err(|e| {
            err(format!(
                "Failed to create output path '{}': {e}",
                output_dir.display()
            ))
        })?;
    }

    let mut output: Vec<(String, (PathBuf, Option<PathBuf>))> = Vec::new();
    for (lane, (f1, f2)) in &undetermined {
        let o1 = output_dir.join(discovery::purged_name(f1));
        let o2 = f2
            .as_ref()
            .map(|f| output_dir.join(discovery::purged_name(f)));
        output.push((lane.clone(), (o1, o2)));
    }

    log.info("Output files:");
    for (lane, (o1, o2)) in &output {
        log.info(&format!("    {lane}: '{}'", o1.display()));
        match o2 {
            Some(o) => log.info(&format!("    {lane}: '{}'", o.display())),
            None => log.info(&format!("    {lane}: 'None'")),
        }
    }

    // ---- assigned ------------------------------------------------------
    let mut assigned_files: Vec<PathBuf> = Vec::new();
    for p in &cli.assigned_path {
        match discovery::explode_path(p, true, &log) {
            Ok(files) => assigned_files.extend(files),
            Err(e) => return Err(logged_err(e)),
        }
    }
    assigned_files.sort();
    assigned_files.dedup();

    // The Python CLI logs the assigned files with the *last undetermined
    // lane* variable (a leftover loop variable, main.py:299). Reproduced
    // here for byte-for-byte log parity.
    let last_undetermined_lane = undetermined.last().expect("non-empty").0.clone();

    let (clean, _skipped) = discovery::clean_assigned(assigned_files, &log);
    let mut assigned: Vec<(String, Vec<PathBuf>)> = Vec::new();
    for f in &clean {
        let Some(name) = f.file_name().and_then(|n| n.to_str()) else {
            continue;
        };
        match discovery::lane_key(name) {
            Some(key) => match assigned.iter_mut().find(|(k, _)| *k == key) {
                Some((_, files)) => files.push(f.clone()),
                None => assigned.push((key, vec![f.clone()])),
            },
            None => log.warning(&format!(
                "File '{}' does not match the expected pattern 'L[0-9]{{3}}'",
                f.display()
            )),
        }
    }

    log.info("Assigned files:");
    for (_key, files) in &assigned {
        for val in files {
            log.info(&format!(
                "    {last_undetermined_lane}: '{}'",
                val.display()
            ));
        }
    }

    Ok(AppConfig {
        log,
        undetermined,
        output,
        output_dir,
        assigned,
        method: cli.method,
        max_items: cli.max_items,
        fpr: cli.fpr,
        threads: cli.threads,
        threading_method: cli.threading_method,
    })
}
