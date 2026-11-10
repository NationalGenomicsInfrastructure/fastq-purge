//! The per-lane pipeline for the SampleSheet input mode (M7.3).
//!
//! Contract: the samplesheet_support branch's per-lane log flow
//! (`main`, L1080-1150), run once per lane in sample-sheet order:
//! `Processing lane` + the three `PosixPath` repr DEBUG lines, the
//! method block, and (exact) the multiprocessing-pool-style lines, then
//! the D6 writer (Writing + Replacing only — never Removed/Kept/saved).
//!
//! Deviations (documented in RUST_REWRITE.md §8):
//!
//! * The branch's sheet-mode approx path is broken for all inputs
//!   (AttributeError on the undetermined tuple in `build_bloom_filter`
//!   and `process_target_file`). Sheet-mode approx here follows the
//!   canonical per-lane M4/D4 contract instead (one bloom per lane,
//!   per-file scan, incremental `new duplicates` counts).
//! * The branch spawns one child process per assigned file (distinct
//!   pids); the port runs rayon tasks and reports a single process id.
//!   Child blocks are flushed in file order (the M5 deterministic-log
//!   rule); the branch prints them in completion order. The branch's
//!   child `Process {pid} started: ... with fork` line and its
//!   per-child memory pair are suppressed in the branch itself (the
//!   pool worker's `multiprocessing` logger stays at the `pool_task`
//!   default 'INFO' — `pool.map` drops the CLI log level), so the port
//!   omits them.
//! * Memory lines come from /proc/meminfo (live values) and the
//!   count-line MB is the set/bloom memory estimate (the branch uses
//!   psutil and `sys.getsizeof`); both are non-deterministic on both
//!   sides and are normalized in the differential.
//! * The branch's uncaught rename error is a traceback; the port emits
//!   the same clean ERROR line as for the caught write error (M1
//!   precedent: clean error where Python traceback'd).

use std::fs;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::{AppConfig, Method};
use crate::fastq::{FastqWriter, FileFastqReader, PairedFastqWriter};
use crate::logging::Logger;
use crate::membership::exact::{hash_name0, ExactSet};
use crate::pipeline::approx;
use crate::pipeline::exact;
use crate::pipeline::{flush_lines, Lines, ScanBlock};

/// Run the per-lane pipeline for a sheet-mode config. Lanes run
/// sequentially in sample-sheet order (the branch's `main` loop); the
/// per-lane assigned-file scan runs on the global pool.
pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    for i in 0..config.undetermined.len() {
        match config.method {
            Method::Exact => exact_lane(&config.log, config, i)?,
            Method::Approx => approx_lane(&config.log, config, i)?,
        }
    }
    config.log.info("Done!");
    Ok(())
}

/// One lane: the branch's `main` loop body for `matching_table[key]`.
fn exact_lane(log: &Logger, config: &AppConfig, i: usize) -> anyhow::Result<()> {
    let key = config.undetermined[i].0.as_str();
    let (u1, u2) = &config.undetermined[i].1;
    let (o1, o2) = &config.output[i].1;
    let files = &config.assigned[i].1;

    let mut lines = Lines::default();
    lane_prologue(&mut lines, key, (u1, u2), (o1, o2), files);

    lines.info("Using exact method");
    lines.append_mem(None);

    lines.info("Building python set from undetermined reads...");
    let (set, set_lines) = match exact::build_set(u1) {
        Ok(x) => x,
        Err((l, e)) => {
            let mut all = lines.into_vec();
            all.extend(l.into_vec());
            flush_lines(log, &all);
            return Err(e);
        }
    };
    lines.append(set_lines);
    lines.append_mem(None);

    // The branch's multiprocessing pool (default mp_pool): the third
    // memory pair is printed at `multiprocessing_process_target_set`
    // start, before the pool INFO lines.
    lines.append_mem(None);
    lines.info("Using multiprocessing pool...");
    lines.info("Purging already assigned reads...");
    lines.debug(format!("Creating {} child processes...", files.len()));

    let pid = std::process::id();
    let results: Vec<Result<(ExactSet, Lines), (Lines, anyhow::Error)>> = files
        .par_iter()
        .map(|f| scan_lane_file(f, &set, pid))
        .collect();

    let mut union = ExactSet::new();
    let mut child_blocks = Vec::new();
    let mut parent_counts = Vec::new();
    for (idx, r) in results.into_iter().enumerate() {
        match r {
            Ok((part, l)) => {
                parent_counts.push((idx, union.merge(part)));
                child_blocks.push(l);
            }
            Err((l, e)) => {
                let mut all = lines.into_vec();
                all.extend(l.into_vec());
                flush_lines(log, &all);
                return Err(e);
            }
        }
    }
    // Child lines in file order (M5 rule; the branch prints them in
    // completion order), then the parent's per-file lines.
    for l in child_blocks {
        lines.append(l);
    }
    for (idx, k) in parent_counts {
        lines.debug(format!(
            "Process finished: {} with {k} already assigned reads",
            exact::file_name(&files[idx])
        ));
    }

    lines.info(format!(
        "Found {} duplicates in the undetermined file.",
        union.len()
    ));
    if union.is_empty() {
        lines.info("Found no duplicates in the undetermined reads. Exiting...");
        flush_lines(log, &lines.into_vec());
        return Ok(());
    }

    finish_write(
        log,
        &mut lines,
        u1,
        u2,
        o1,
        o2,
        &union,
        config.keep_original,
    )
}

/// Sheet-mode approx: the canonical per-lane M4/D4 contract (the
/// branch's approx path is broken for all inputs — see module docs).
fn approx_lane(log: &Logger, config: &AppConfig, i: usize) -> anyhow::Result<()> {
    let key = config.undetermined[i].0.as_str();
    let (u1, u2) = &config.undetermined[i].1;
    let (o1, o2) = &config.output[i].1;
    let files = &config.assigned[i].1;

    let mut lines = Lines::default();
    lane_prologue(&mut lines, key, (u1, u2), (o1, o2), files);

    lines.info("Using bloom filter method");
    let (bloom, bloom_lines) = match approx::build_bloom(config.max_items, config.fpr, u1) {
        Ok(x) => x,
        Err((l, e)) => {
            let mut all = lines.into_vec();
            all.extend(l.into_vec());
            flush_lines(log, &all);
            return Err(e);
        }
    };
    lines.append(bloom_lines);

    lines.info("Scanning assigned files for already-assigned reads...");
    let results: Vec<ScanBlock> = files
        .par_iter()
        .map(|f| approx::scan_file(f, std::slice::from_ref(&bloom)))
        .collect();

    let mut union = ExactSet::new();
    let mut blocks = Vec::new();
    for (idx, r) in results.into_iter().enumerate() {
        match r {
            // One bloom => one part: the lane's own partial set.
            Ok((parts, l)) => {
                let k = union.merge(parts.into_iter().next().expect("one part per bloom"));
                blocks.push((idx, l, k));
            }
            Err((l, e)) => {
                let mut all = lines.into_vec();
                all.extend(l.into_vec());
                flush_lines(log, &all);
                return Err(e);
            }
        }
    }
    for (idx, l, k) in blocks {
        lines.append(l);
        lines.debug(format!(
            "'{}': {k} new duplicates",
            exact::file_name(&files[idx])
        ));
    }

    lines.info(format!(
        "Found {} duplicates in the undetermined file.",
        union.len()
    ));
    if union.is_empty() {
        lines.info("Found no duplicates in the undetermined reads. Exiting...");
        flush_lines(log, &lines.into_vec());
        return Ok(());
    }

    finish_write(
        log,
        &mut lines,
        u1,
        u2,
        o1,
        o2,
        &union,
        config.keep_original,
    )
}

/// The branch's per-lane prologue: `Processing lane` + the three file
/// repr DEBUG lines (Python list/tuple reprs of PosixPath values).
fn lane_prologue(
    lines: &mut Lines,
    key: &str,
    und: (&PathBuf, &Option<PathBuf>),
    output: (&PathBuf, &Option<PathBuf>),
    files: &[PathBuf],
) {
    lines.info(format!("Processing lane '{key}'"));
    lines.debug(format!(
        "    Undetermined files: [({})]",
        pair_repr(und.0, und.1)
    ));
    lines.debug(format!(
        "    Output files: ({})",
        pair_repr(output.0, output.1)
    ));
    lines.debug(format!(
        "    Assigned files: [{}]",
        files.iter().map(|p| pp(p)).collect::<Vec<_>>().join(", ")
    ));
}

/// The shared write tail: the D6 writer, then flush + (on write failure)
/// the branch's ERROR line and a suppressed error (already reported).
// Eight parameters: the lane's four paths + the writer's inputs; a
// wrapper struct would only move the arguments into one field.
#[allow(clippy::too_many_arguments)]
fn finish_write(
    log: &Logger,
    lines: &mut Lines,
    u1: &Path,
    u2: &Option<PathBuf>,
    o1: &Path,
    o2: &Option<PathBuf>,
    dups: &ExactSet,
    keep: bool,
) -> anyhow::Result<()> {
    match write_sheet(u1, u2, o1, o2, dups, keep) {
        Ok(wl) => {
            lines.append(wl);
            flush_lines(log, &lines.to_vec());
            Ok(())
        }
        Err((wl, msg)) => {
            let mut all = lines.to_vec();
            all.extend(wl.into_vec());
            flush_lines(log, &all);
            log.error(&msg);
            // Reported through the logger: main must not print it again.
            Err(anyhow::anyhow!(crate::pipeline::Logged))
        }
    }
}

/// The branch's `write_purged_fastq` (D6): the Writing line(s), the
/// filtered write, and — unless `keep` — the Replacing line (trailing
/// space, branch quirk) plus the in-place rename that destroys the
/// original. Unlike the legacy writer it never emits Removed/Kept/saved.
///
/// On failure returns the buffered lines plus the branch's ERROR
/// message.
fn write_sheet(
    u1: &Path,
    u2: &Option<PathBuf>,
    o1: &Path,
    o2: &Option<PathBuf>,
    dups: &ExactSet,
    keep: bool,
) -> Result<Lines, (Lines, String)> {
    let mut lines = Lines::default();
    match (u2, o2) {
        (Some(u2), Some(o2)) => {
            lines.info(format!(
                "Writing purged fastq files '{}' and '{}'",
                exact::file_name(o1),
                exact::file_name(o2)
            ));
            let w = (|| -> anyhow::Result<()> {
                let mut reader1 = FileFastqReader::open(u1)?;
                let mut reader2 = FileFastqReader::open(u2)?;
                let mut writer = PairedFastqWriter::create(o1, o2)?;
                while let Some(rec) = reader1.next_record() {
                    let rec1 = rec?;
                    // zip semantics: stop at the shorter of R1/R2
                    let rec2 = match reader2.next_record() {
                        Some(r) => r?,
                        None => break,
                    };
                    if !dups.contains_name(rec1.name0()) {
                        writer.write_pair(&rec1, &rec2)?;
                    }
                }
                writer.finish()?;
                Ok(())
            })();
            if let Err(e) = w {
                return Err((
                    lines,
                    format!(
                        "Error writing purged fastq files '{}' and '{}': {e}",
                        exact::file_name(o1),
                        exact::file_name(o2)
                    ),
                ));
            }
            if !keep {
                lines.info(format!(
                    "Replacing undetermined fastq files '{}' and '{}' with purged fastq files ",
                    exact::file_name(u1),
                    exact::file_name(u2)
                ));
                if let Err(e) = fs::rename(o1, u1).and_then(|_| fs::rename(o2, u2)) {
                    return Err((
                        lines,
                        format!(
                            "Error writing purged fastq files '{}' and '{}': {e}",
                            exact::file_name(o1),
                            exact::file_name(o2)
                        ),
                    ));
                }
            }
        }
        _ => {
            // Single-read lane (u2/o2 None; the output tuple mirrors the
            // input tuple, so a mixed shape cannot occur).
            lines.info(format!(
                "Writing purged fastq file '{}'",
                exact::file_name(o1)
            ));
            let w = (|| -> anyhow::Result<()> {
                let mut reader = FileFastqReader::open(u1)?;
                let mut writer = FastqWriter::create(o1)?;
                while let Some(rec) = reader.next_record() {
                    let rec = rec?;
                    if !dups.contains_name(rec.name0()) {
                        writer.write_record(&rec.name, &rec.seq, &rec.qual)?;
                    }
                }
                writer.finish()?;
                Ok(())
            })();
            if let Err(e) = w {
                return Err((
                    lines,
                    format!(
                        "Error writing purged fastq file '{}': {e}",
                        exact::file_name(o1)
                    ),
                ));
            }
            if !keep {
                lines.info(format!(
                    "Replacing undetermined fastq file '{}' ",
                    exact::file_name(u1)
                ));
                if let Err(e) = fs::rename(o1, u1) {
                    return Err((
                        lines,
                        format!(
                            "Error writing purged fastq file '{}': {e}",
                            exact::file_name(o1)
                        ),
                    ));
                }
            }
        }
    }
    Ok(lines)
}

/// The branch's per-child `pool_task` line for one assigned file:
/// `Process {pid} ('{name}') finished yielding {k} already assigned
/// reads`. `k` is the absolute number of distinct names in this file
/// that are in the lane set (the branch's per-child result).
///
/// The branch's child `Process {pid} started: ... with fork` DEBUG line
/// and its per-child memory pair go through a `multiprocessing` logger
/// left at the `pool_task` default level 'INFO' (`pool.map` does not
/// propagate the CLI log level), so the branch suppresses them even for
/// DEBUG runs; they are not emitted here (observable-parity port).
fn scan_lane_file(
    f: &Path,
    set: &ExactSet,
    pid: u32,
) -> Result<(ExactSet, Lines), (Lines, anyhow::Error)> {
    let name = exact::file_name(f);
    let mut lines = Lines::default();

    let mut part = ExactSet::new();
    let result = (|| -> anyhow::Result<()> {
        let mut reader = FileFastqReader::open(f)?;
        while let Some(rec) = reader.next_record() {
            let rec = rec?;
            let h = hash_name0(rec.name0());
            if set.contains_hash(h) {
                part.insert_hash(h);
            }
        }
        Ok(())
    })();
    if let Err(e) = result {
        return Err((lines, e));
    }

    let k = part.len();
    lines.debug(format!(
        "Process {pid} ('{name}') finished yielding {k} already assigned reads"
    ));
    Ok((part, lines))
}

/// A Python `PosixPath` repr (the branch's paths are PosixPath on
/// Linux).
fn pp(p: &Path) -> String {
    format!("PosixPath('{}')", p.display())
}

/// A Python tuple repr of an (R1, R2 or None) pair.
fn pair_repr(p1: &Path, p2: &Option<PathBuf>) -> String {
    match p2 {
        Some(p2) => format!("{}, {}", pp(p1), pp(p2)),
        None => format!("{}, None", pp(p1)),
    }
}

/// The branch's `memory_usage` pair, from /proc/meminfo (live values —
/// non-deterministic on both sides; the differential normalizes these
/// lines). `pid` gives the `PID {pid} - ` prefix (the branch's
/// per-child lines).
pub(crate) fn mem_lines(pid: Option<u32>) -> [String; 2] {
    let m = read_meminfo();
    let mem_used = m.total.saturating_sub(m.available);
    // `mem_used` is 0 whenever `total` is 0 (saturating sub), so `max(1)`
    // keeps the 0% result for a bogus /proc/meminfo.
    let mem_pct = mem_used * 100 / m.total.max(1);
    let swap_used = m.swap_total.saturating_sub(m.swap_free);
    let swap_pct = swap_used * 100 / m.swap_total.max(1);
    let prefix = match pid {
        Some(pid) => format!("PID {pid} - "),
        None => String::new(),
    };
    [
        format!(
            "{prefix}Memory used: {:.2} GB ({mem_pct}%)",
            mem_used as f64 / 1024.0 / 1024.0
        ),
        format!(
            "{prefix}Swap used: {:.2} GB ({swap_pct}%)",
            swap_used as f64 / 1024.0 / 1024.0
        ),
    ]
}

#[derive(Default)]
struct Meminfo {
    total: u64,     // kB
    available: u64, // kB
    swap_total: u64,
    swap_free: u64,
}

fn read_meminfo() -> Meminfo {
    let mut m = Meminfo::default();
    if let Ok(s) = fs::read_to_string("/proc/meminfo") {
        for line in s.lines() {
            if let Some((k, v)) = line.split_once(':') {
                let kb = v
                    .trim()
                    .trim_end_matches("kB")
                    .trim()
                    .parse::<u64>()
                    .unwrap_or(0);
                match k.trim() {
                    "MemTotal" => m.total = kb,
                    "MemAvailable" => m.available = kb,
                    "SwapTotal" => m.swap_total = kb,
                    "SwapFree" => m.swap_free = kb,
                    _ => {}
                }
            }
        }
    }
    m
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;

    /// A lane's files in a temp dir: undetermined pair, assigned pair
    /// (one name shared with the undetermined set), and the expected
    /// purged output.
    fn setup(tag: &str) -> (PathBuf, PathBuf, PathBuf, PathBuf, PathBuf, PathBuf) {
        let dir = std::env::temp_dir().join(format!("fqp_m73_{tag}_{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        // Plain (non-gz) names: the reader/writer dispatch on extension.
        let und1 = dir.join("Undetermined_S1_L001_R1_001.fastq");
        let und2 = dir.join("Undetermined_S1_L001_R2_001.fastq");
        let out1 = dir.join("Undetermined_S1_L001_R1_001.purged.fastq");
        let out2 = dir.join("Undetermined_S1_L001_R2_001.purged.fastq");
        let ass = dir.join("S1_L001_R1_001.fastq");
        write_fastq(&und1, "@u1a 1\nAAAA\n+\n*@@@\n@u1b 1\nTTTT\n+\n*@@@\n");
        write_fastq(&und2, "@u1a/2 1\nAAAA\n+\n*@@@\n@u1b/2 1\nTTTT\n+\n*@@@\n");
        // a1a is a duplicate of u1a; a1b is not.
        write_fastq(&ass, "@a1a 1\nCCCC\n+\n*@@@\n@a1b 1\nGGGG\n+\n*@@@\n");
        (dir, und1, und2, out1, out2, ass)
    }

    fn write_fastq(p: &Path, body: &str) {
        let mut f = fs::File::create(p).unwrap();
        f.write_all(body.as_bytes()).unwrap();
    }

    fn dups() -> ExactSet {
        // name0 of "@u1a 1" is "u1a" (Record.name excludes the leading '@').
        let mut s = ExactSet::new();
        s.insert_hash(hash_name0(b"u1a"));
        s
    }

    #[test]
    fn write_sheet_replace_and_keep() {
        // !keep: the .purged output replaces the original in place.
        let (_dir, und1, und2, out1, out2, _ass) = setup("replace");
        let d = dups();
        let u2 = Some(und2.clone());
        let o2 = Some(out2.clone());
        let res = write_sheet(&und1, &u2, &out1, &o2, &d, false);
        assert!(res.is_ok());
        // The rename moved the purged files onto the original paths.
        assert!(!out1.exists());
        assert!(!out2.exists());
        let body = fs::read_to_string(&und1).unwrap();
        assert!(body.contains("@u1b 1"));
        assert!(!body.contains("@u1a 1"));

        // keep: the purged files stay next to the originals.
        let (_dir, und1, und2, out1, out2, _ass) = setup("keep");
        let d = dups();
        let u2 = Some(und2.clone());
        let o2 = Some(out2.clone());
        let res = write_sheet(&und1, &u2, &out1, &o2, &d, true);
        assert!(res.is_ok());
        assert!(out1.exists());
        assert!(und1.exists());
        let body = fs::read_to_string(&out1).unwrap();
        assert!(body.contains("@u1b 1"));
        assert!(!body.contains("@u1a 1"));
    }
}
