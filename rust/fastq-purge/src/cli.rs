//! CLI definitions and validation entry point (milestone M1).
//!
//! Flags, defaults and choices mirror `fastq_purge.main.parse_args()`
//! (Python, argparse) exactly; discovery/validation behaviour mirrors
//! `validate_args()`. See RUST_REWRITE.md §7 for the parity checklist.

use std::path::{Path, PathBuf};

use clap::Parser;

use crate::discovery::{self, UndeterminedLane};
use crate::logging::{Level, Logger};
use crate::samplesheet;
use crate::sheetdiscovery;

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
    /// Path to the flowcell directory. It is used to extract the lane
    /// information from the SampleSheet.csv file, and to identify the
    /// undetermined and assigned files.
    #[arg(long)]
    pub flowcell_path: Option<PathBuf>,

    /// Path to the target fastq file(s) to purge. It can be gzipped or not.
    /// It can be a single file or a directory. If a directory is provided,
    /// all files in the directory that match the pattern '*.fq*' or
    /// '*.fastq*' will be used as target files. The search is non-recursive.
    #[arg(long)]
    pub undetermined_path: Option<PathBuf>,

    /// Path to the output fastq file
    #[arg(long)]
    pub output_path: Option<String>,

    /// Path to the bloom filter sources file(s). They can be gzipped or
    /// not. It can be a single file or multiple files separated by spaces,
    /// or a directory. If a directory is provided, all files in the
    /// directory that match the pattern '*.fq*' or '*.fastq*' will be used
    /// as bloom filter sources. The search is recursive.
    #[arg(long, num_args = 1..)]
    pub assigned_path: Vec<PathBuf>,

    /// Path to the sample sheet file. It is used to extract the sample and
    /// lane information, and to match the undetermined files with the
    /// assigned files.
    #[arg(long)]
    pub sample_sheet: Option<PathBuf>,

    /// Keep the original undetermined fastq files after purging. If this
    /// option is not set, the original files will be removed after purging,
    /// while if set, the purged files will be saved with the suffix
    /// '.purged' added to the original file basename.
    #[arg(long)]
    pub keep_original: bool,

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

    /// Purge the assigned files instead of the undetermined ones.
    ///
    /// Canonical direction (default): undetermined reads whose name
    /// appears in any assigned file are removed (RUST_REWRITE.md D1).
    /// This flag inverts that for anyone relying on the old approx
    /// direction (D3).
    #[arg(long, default_value_t = false)]
    pub purge_assigned: bool,
}

/// Validated, ready-to-run configuration (mirrors the Python
/// `args` namespace after `validate_args`).
///
/// The pipeline fields are consumed from M3/M4 onward; M1's stub
/// pipeline reads none of them.
/// The input mode selected by the CLI (M7, decision D8):
///
/// * `Legacy` — no `--flowcell-path` and no `--sample-sheet`: the M1-M5
///   directory-discovery behaviour, byte-identical.
/// * `Sheet` — the SampleSheet alternate input mode (M7): per-lane
///   pipeline; the undetermined originals are replaced in place unless
///   `--keep-original` is set (decision D6).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunMode {
    Legacy,
    Sheet,
}

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
    /// lane key -> assigned source files, in discovery order. Sheet mode:
    /// the per-project inner lists flattened in sheet order (empty inner
    /// lists included, as in the branch).
    pub assigned: Vec<(String, Vec<PathBuf>)>,
    pub method: Method,
    pub max_items: u64,
    pub fpr: f64,
    pub threads: u64,
    pub threading_method: ThreadingMethod,
    /// D3: invert the canonical purge direction (purge assigned instead of
    /// undetermined).
    pub purge_assigned: bool,
    /// The selected input mode (M7, D8).
    pub mode: RunMode,
    /// Sheet mode: keep the original undetermined files after purging
    /// (D6; ignored in legacy mode).
    pub keep_original: bool,
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

    // ---- M7, D8: input-mode dispatch ------------------------------------
    // * `--flowcell-path` (+ optional overrides): flowcell mode.
    // * `--sample-sheet` without `--flowcell-path`: manual sheet mode
    //   (D5: the branch crashes here with AttributeError on
    //   `None.is_dir()`; we implement the intended else-branch).
    // * Neither: the legacy M1-M5 directory-discovery mode.
    let fc = cli.flowcell_path.clone();
    let sheet_arg = cli.sample_sheet.clone();
    let und_arg = cli.undetermined_path.clone();
    let assigned_arg = cli.assigned_path.clone();

    if let Some(fc) = fc {
        if !fc.is_dir() {
            return Err(logged_err(format!(
                "Flowcell path '{}' is not a valid directory! Please provide a valid flowcell path.",
                fc.display()
            )));
        }
        log.debug("Flowcell path provided! Validating paths...");
        let sheet = match &sheet_arg {
            Some(s) => s.clone(),
            None => find_flowcell_sheet(&fc).ok_or_else(|| {
                // The branch indexes [-1] of an empty list (IndexError);
                // we fail cleanly (M1 precedent).
                err(format!(
                    "No sample sheet found in flowcell '{}'. Please provide a valid sample sheet file.",
                    fc.display()
                ))
            })?,
        };
        let assigned = if assigned_arg.is_empty() {
            vec![fc.clone()]
        } else {
            assigned_arg
        };
        let und = und_arg.unwrap_or(fc);
        return validate_sheet_mode(cli, log, sheet, Some(und), assigned);
    }
    if let Some(sheet) = sheet_arg {
        log.info("No flowcell path provided, using manually provided arguments");
        return validate_sheet_mode(cli, log, sheet, und_arg, assigned_arg);
    }
    validate_legacy(cli, log)
}

/// Sheet-mode validation: the branch's `validate_args` common checks +
/// `parse_sample_sheet` + `process_args` discovery (sheetdiscovery).
fn validate_sheet_mode(
    cli: Cli,
    log: Logger,
    sheet: PathBuf,
    und: Option<PathBuf>,
    assigned: Vec<PathBuf>,
) -> Result<AppConfig, (String, i32, bool)> {
    if !sheet.is_file() {
        return Err(logged_err(
            "No sample sheet provided or found! Please provide a valid sample sheet file."
                .to_string(),
        ));
    }
    let und = match und {
        Some(p) if p.exists() => p,
        _ => {
            return Err(logged_err(
                "No undetermined path provided or found! Please provide a valid undetermined fastq file."
                    .to_string(),
            ))
        }
    };
    if assigned.is_empty() || !assigned.iter().all(|p| p.exists()) {
        return Err(logged_err(
            "No assigned path provided or found! Please provide assigned fastq file(s) or directory."
                .to_string(),
        ));
    }
    log.debug("Success! All paths are valid.");

    if cli.purge_assigned {
        return Err(err(
            "--purge-assigned is only supported in directory-discovery mode (it is not part of the SampleSheet input mode).".to_string(),
        ));
    }

    let data = samplesheet::parse(&sheet).map_err(|e| (format!("{e}"), 1, false))?;
    let sheet_plan = sheetdiscovery::plan(
        &und,
        &assigned,
        cli.output_path.as_deref().map(Path::new),
        &data,
        &log,
    )
    .map_err(|(msg, logged)| (msg, 1, logged))?;

    let undetermined: Vec<(String, UndeterminedLane)> = sheet_plan
        .lanes
        .iter()
        .map(|l| (l.key.clone(), l.undetermined.clone()))
        .collect();
    let output: Vec<(String, (PathBuf, Option<PathBuf>))> = sheet_plan
        .lanes
        .iter()
        .map(|l| (l.key.clone(), l.output.clone()))
        .collect();
    let assigned: Vec<(String, Vec<PathBuf>)> = sheet_plan
        .lanes
        .iter()
        .map(|l| {
            (
                l.key.clone(),
                l.assigned.iter().flatten().cloned().collect(),
            )
        })
        .collect();

    Ok(AppConfig {
        log,
        undetermined,
        output,
        output_dir: sheet_plan.output_dir,
        assigned,
        method: cli.method,
        max_items: cli.max_items,
        fpr: cli.fpr,
        threads: cli.threads,
        threading_method: cli.threading_method,
        purge_assigned: false,
        mode: RunMode::Sheet,
        keep_original: cli.keep_original,
    })
}

/// The last flowcell entry (directory-entry order) whose name matches the
/// branch regex `^SampleSheet.(csv|txt)$` — unescaped dots, so the two dot
/// positions are any single character.
fn find_flowcell_sheet(fc: &Path) -> Option<PathBuf> {
    let mut last = None;
    for entry in std::fs::read_dir(fc).ok()?.flatten() {
        let file_name = entry.file_name();
        let name = file_name.to_string_lossy();
        if is_sample_sheet_name(&name) {
            last = Some(entry.path());
        }
    }
    last
}

fn is_sample_sheet_name(name: &str) -> bool {
    let b = name
        .strip_prefix("SampleSheet")
        .unwrap_or_default()
        .as_bytes();
    b.len() == 4 && (b[1..] == *b"csv" || b[1..] == *b"txt")
}

/// Legacy (M1-M5) directory-discovery validation, byte-identical to the
/// pre-M7 behaviour.
fn validate_legacy(cli: Cli, log: Logger) -> Result<AppConfig, (String, i32, bool)> {
    // Usage errors (rc 2). The path flags are optional at the parser level
    // since M7, so the required-argument check is manual here, reproducing
    // clap's error text captured pre-M7.
    let mut missing = Vec::new();
    if cli.undetermined_path.is_none() {
        missing.push("  --undetermined-path <UNDETERMINED_PATH>".to_string());
    }
    if cli.assigned_path.is_empty() {
        missing.push("  --assigned-path <ASSIGNED_PATH>...".to_string());
    }
    if !missing.is_empty() {
        return Err((
            format!(
                "error: the following required arguments were not provided:\n{}\n\nUsage: fastq-purge-rs [options]\n\nFor more information, try '--help'.",
                missing.join("\n")
            ),
            2,
            false,
        ));
    }
    let und_path = cli.undetermined_path.clone().unwrap();

    // ---- undetermined -------------------------------------------------
    let un_files = discovery::explode_path(&und_path, false, &log).map_err(logged_err)?;

    let undetermined = discovery::create_paired_dict(un_files).map_err(err)?;
    if undetermined.is_empty() {
        // The Python CLI raises IndexError here (empty dict -> keys()[0]);
        // we fail cleanly instead.
        return Err(err(format!(
            "No undetermined files found in '{}' (pattern '*.fq*' or '*.fastq*')",
            und_path.display()
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
        purge_assigned: cli.purge_assigned,
        mode: RunMode::Legacy,
        keep_original: false,
    })
}
