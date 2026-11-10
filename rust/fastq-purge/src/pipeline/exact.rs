//! The exact pipeline (milestone M3; M5: rayon), canonical semantics per
//! RUST_REWRITE.md §3/§5 (decisions D1–D3, accepted 2026-11-06):
//!
//! - D1: the undetermined files are purged — an undetermined read is
//!   removed when its name0 appears in any assigned file (the old Python
//!   exact direction). `--purge-assigned` (D3) inverts this.
//! - D2: the comparison key is name0 (header up to the first space) in
//!   both methods.
//! - Single scan: the assigned files are streamed once and each name0 is
//!   tested against all per-lane undetermined sets (the Python code
//!   re-scans the assigned files once per undetermined lane).
//!
//! Observable outputs on the same inputs match the Python exact method
//! (same purged files, same "Found N duplicates" counts); the log
//! sequence is the new canonical contract below (the Python sequence is
//! per-lane: build → scan-all-assigned → write, repeated K times).
//!
//! M5 (rayon): the passes run on the global pool sized by `--threads`
//! (main.rs). Deterministic-log rule: workers never print — each returns
//! the lines the sequential run would have emitted for its item (a
//! `Lines` block), and the orchestrator flushes blocks in canonical
//! lane/file order. The log is therefore byte-identical to the M3
//! sequential run for any `--threads` (and any failure: the block of the
//! failing item is flushed up to the point of failure, then the error is
//! returned). Pass 2 merges per-file partial dups sets in file order, so
//! the per-file `new duplicates` count equals the sequential
//! stream-order count exactly (see `ExactSet::merge`).
//!
//! Log contract (INFO level):
//! ```text
//! Using exact method
//! Building membership set from undetermined reads...
//!   (DEBUG) Adding '{und R1 file name}' to set            [per lane]
//!   Number of undetermined reads: {N} ({MB:.2f} MB)       [per lane]
//! Scanning assigned files for already-assigned reads...
//!   Processing target file '{assigned file name}'         [per file]
//!   (DEBUG) '{assigned file name}': {k} new duplicates    [per file]
//! Found {K} duplicates in the undetermined file.          [per lane]
//! [lane with K > 0:]
//!   Writing purged fastq file '{out1 name}'
//!   | Writing purged fastq files '{out1 name}' and '{out2 name}'
//!       Removed {R} reads
//!       Kept {W} reads
//!       Purged fastq file saved to '{out path}'
//! [lane with K == 0:]
//!   Found no duplicates in the undetermined reads. Exiting...
//! Done!
//! ```
//! (`--purge-assigned` replaces the scan/write block with per-assigned-
//! file lines; see `purge_assigned`.)

use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::AppConfig;
use crate::discovery::purged_name;
use crate::fastq::{FastqWriter, FileFastqReader, PairedFastqWriter};
use crate::membership::exact::{hash_name0, ExactSet};
use crate::pipeline::{flush_lines, Block, Lines, ScanBlock};

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    config.log.info("Using exact method");

    // ---- pass 1: membership set per undetermined lane ----------------
    config
        .log
        .info("Building membership set from undetermined reads...");
    let built: Vec<Result<(ExactSet, Lines), (Lines, anyhow::Error)>> = config
        .undetermined
        .par_iter()
        .map(|(_lane, (r1, _r2))| build_set(r1))
        .collect();
    let mut sets: Vec<ExactSet> = Vec::with_capacity(built.len());
    for b in built {
        match b {
            Ok((set, lines)) => {
                flush_lines(&config.log, &lines.v);
                sets.push(set);
            }
            Err((lines, e)) => {
                flush_lines(&config.log, &lines.v);
                return Err(e);
            }
        }
    }

    if config.purge_assigned {
        return purge_assigned(config, &sets);
    }

    // ---- pass 2: one scan of the assigned files, all lane sets --------
    config
        .log
        .info("Scanning assigned files for already-assigned reads...");
    let k = sets.len();
    let files: Vec<&Path> = config
        .assigned
        .iter()
        .flat_map(|(_lane, fs)| fs.iter().map(PathBuf::as_path))
        .collect();
    let scanned: Vec<ScanBlock> = files.par_iter().map(|f| scan_file(f, &sets)).collect();
    let mut dups: Vec<ExactSet> = (0..k).map(|_| ExactSet::new()).collect();
    for (i, s) in scanned.into_iter().enumerate() {
        let (mut parts, lines) = match s {
            Ok(v) => v,
            Err((lines, e)) => {
                flush_lines(&config.log, &lines.v);
                return Err(e);
            }
        };
        // canonical position: the file's Processing line, then its
        // merge count (the count needs the sequential merge state)
        flush_lines(&config.log, &lines.v);
        let mut new_here = 0usize;
        for l in 0..k {
            new_here += dups[l].merge(std::mem::take(&mut parts[l]));
        }
        config.log.debug(&format!(
            "'{}': {new_here} new duplicates",
            file_name(files[i])
        ));
    }

    // ---- pass 3: write the purged undetermined files -------------------
    let n_lanes = dups.len();
    let writes: Vec<Option<Block>> = (0..n_lanes)
        .into_par_iter()
        .map(|l| {
            if dups[l].is_empty() {
                return None;
            }
            let (_lane, (r1, r2)) = &config.undetermined[l];
            let (_lane, (out1, out2)) = &config.output[l];
            Some(write_purged_undetermined(
                r1,
                r2.as_deref(),
                out1,
                out2.as_deref(),
                &dups[l],
                false,
            ))
        })
        .collect();
    for (l, w) in writes.into_iter().enumerate() {
        let d = dups[l].len();
        config
            .log
            .info(&format!("Found {d} duplicates in the undetermined file."));
        if d == 0 {
            config
                .log
                .info("Found no duplicates in the undetermined reads. Exiting...");
            continue;
        }
        match w.expect("worker dispatched for every lane with dups") {
            Ok(lines) => flush_lines(&config.log, &lines.v),
            Err((lines, e)) => {
                flush_lines(&config.log, &lines.v);
                return Err(e);
            }
        }
    }

    config.log.info("Done!");
    Ok(())
}

/// Pass-1 worker: stream one undetermined R1 file into an `ExactSet`.
/// Returns the set plus the log lines in sequential emission order
/// (`Adding ... to set` first, then the count line).
pub(crate) fn build_set(r1: &Path) -> Result<(ExactSet, Lines), (Lines, anyhow::Error)> {
    let mut lines = Lines::default();
    lines.debug(format!("Adding '{}' to set", file_name(r1)));
    let mut set = ExactSet::new();
    let mut reader = FileFastqReader::open(r1).map_err(|e| (lines.clone(), e.into()))?;
    while let Some(name0) = reader.next_name0() {
        set.insert_name(name0.map_err(|e| (lines.clone(), e.into()))?);
    }
    let n = set.len();
    let mb = set.estimated_bytes() as f64 / 1024.0 / 1024.0;
    lines.info(format!("Number of undetermined reads: {n} ({mb:.2} MB)"));
    Ok((set, lines))
}

/// Pass-2 worker: stream one assigned file and collect, per undetermined
/// lane, the partial set of lane-member name0 hashes seen in this file
/// (insertion is idempotent — the per-file new-duplicate count is
/// computed by the orchestrator at merge time, see `run`).
pub(crate) fn scan_file(f: &Path, sets: &[ExactSet]) -> ScanBlock {
    let k = sets.len();
    let mut lines = Lines::default();
    lines.info(format!("Processing target file '{}'", file_name(f)));
    let mut parts: Vec<ExactSet> = (0..k).map(|_| ExactSet::new()).collect();
    let mut reader = FileFastqReader::open(f).map_err(|e| (lines.clone(), e.into()))?;
    while let Some(name0) = reader.next_name0() {
        let h = hash_name0(name0.map_err(|e| (lines.clone(), e.into()))?);
        for l in 0..k {
            if sets[l].contains_hash(h) {
                parts[l].insert_hash(h);
            }
        }
    }
    Ok((parts, lines))
}

/// D3: purge the ASSIGNED files instead (assigned reads whose name0
/// appears in the lane's undetermined set are removed). Assigned files
/// are matched to undetermined lanes by lane key; an assigned lane with
/// no undetermined lane is skipped with a warning. M5: the per-file
/// workers run in parallel; the blocks are flushed in lane/file order.
fn purge_assigned(config: &AppConfig, sets: &[ExactSet]) -> anyhow::Result<()> {
    config.log.info("Purging assigned files...");

    // warn once per assigned lane that has no undetermined counterpart
    let mut warned: Vec<String> = Vec::new();
    for (key, _files) in &config.assigned {
        if !config.undetermined.iter().any(|(l, _)| l == key) && !warned.contains(key) {
            warned.push(key.clone());
            config.log.warning(&format!(
                "Lane {key}: no undetermined files; its assigned files are not purged"
            ));
        }
    }

    // per-lane work in canonical (undetermined) lane order
    let lane_files: Vec<Option<Vec<&Path>>> = config
        .undetermined
        .iter()
        .map(|(lane, _und)| {
            config
                .assigned
                .iter()
                .find(|(k, _)| k == lane)
                .map(|(_key, files)| files.iter().map(PathBuf::as_path).collect())
        })
        .collect();
    let mut results: Vec<Option<Vec<Block>>> = lane_files
        .par_iter()
        .enumerate()
        .map(|(l, lf)| {
            lf.as_ref().map(|files| {
                files
                    .par_iter()
                    .map(|f| purge_file(f, &config.output_dir, &sets[l]))
                    .collect()
            })
        })
        .collect();
    for (l, (lane, _und)) in config.undetermined.iter().enumerate() {
        match results[l].take() {
            None => config
                .log
                .debug(&format!("lane {lane}: no assigned files to purge")),
            Some(blocks) => {
                for b in blocks {
                    match b {
                        Ok(lines) => flush_lines(&config.log, &lines.v),
                        Err((lines, e)) => {
                            flush_lines(&config.log, &lines.v);
                            return Err(e);
                        }
                    }
                }
            }
        }
    }

    config.log.info("Done!");
    Ok(())
}

/// `--purge-assigned` worker: one assigned file. The output path is
/// `--output-path / purged_name(file)`.
fn purge_file(f: &Path, out_dir: &Path, set: &ExactSet) -> Block {
    let out = out_dir.join(purged_name(f));
    let mut lines = Lines::default();
    lines.info(format!("Processing target file '{}'", file_name(f)));
    lines.info(format!("    Output file: '{}'", out.display()));
    let mut reader = FileFastqReader::open(f).map_err(|e| (lines.clone(), e.into()))?;
    let mut writer = FastqWriter::create(&out).map_err(|e| (lines.clone(), e.into()))?;
    let (mut kept, mut removed) = (0u64, 0u64);
    while let Some(rec) = reader.next_record() {
        let rec = rec.map_err(|e| (lines.clone(), e.into()))?;
        if set.contains_name(rec.name0()) {
            removed += 1;
        } else {
            writer
                .write_record(&rec.name, &rec.seq, &rec.qual)
                .map_err(|e| (lines.clone(), e.into()))?;
            kept += 1;
        }
    }
    writer.finish().map_err(|e| (lines.clone(), e.into()))?;
    lines.info(format!("    Removed {removed} reads"));
    lines.info(format!("    Kept {kept} reads"));
    lines.info(format!(
        "    Purged fastq file saved to '{}'",
        out.display()
    ));
    Ok(lines)
}

/// Write the purged undetermined file(s) of one lane (the default D1
/// direction). Paired inputs are zipped: iteration stops at the shorter
/// of R1/R2 (Python `dnaio` paired iteration is a zip). `progress`
/// collects the per-1M `    Processed {n} reads` line (the approx
/// method's contract; the exact method does not print it — Python's
/// exact path doesn't either). M5: returns the log lines instead of
/// printing them (the orchestrator flushes in lane order). Shared by the
/// exact and approx pipelines (M4: the dups set is exact in both —
/// decision D4).
pub(crate) fn write_purged_undetermined(
    r1: &Path,
    r2: Option<&Path>,
    out1: &Path,
    out2: Option<&Path>,
    dups: &ExactSet,
    progress: bool,
) -> Block {
    let mut lines = Lines::default();
    let (mut kept, mut removed) = (0u64, 0u64);
    let mut seen = 0u64;

    match (r2, out2) {
        (Some(r2), Some(out2)) => {
            lines.info(format!(
                "Writing purged fastq files '{}' and '{}'",
                file_name(out1),
                file_name(out2)
            ));
            let mut reader1 = FileFastqReader::open(r1).map_err(|e| (lines.clone(), e.into()))?;
            let mut reader2 = FileFastqReader::open(r2).map_err(|e| (lines.clone(), e.into()))?;
            let mut writer =
                PairedFastqWriter::create(out1, out2).map_err(|e| (lines.clone(), e.into()))?;
            while let Some(rec1) = reader1.next_record() {
                let rec1 = rec1.map_err(|e| (lines.clone(), e.into()))?;
                // zip semantics: stop at the shorter of R1/R2
                let rec2 = match reader2.next_record() {
                    Some(r) => r.map_err(|e| (lines.clone(), e.into()))?,
                    None => break,
                };
                seen += 1;
                if progress && seen.is_multiple_of(1_000_000) {
                    lines.info(format!("    Processed {seen} reads"));
                }
                if dups.contains_name(rec1.name0()) {
                    removed += 1;
                } else {
                    writer
                        .write_pair(&rec1, &rec2)
                        .map_err(|e| (lines.clone(), e.into()))?;
                    kept += 1;
                }
            }
            writer.finish().map_err(|e| (lines.clone(), e.into()))?;
        }
        (None, None) => {
            lines.info(format!("Writing purged fastq file '{}'", file_name(out1)));
            let mut reader = FileFastqReader::open(r1).map_err(|e| (lines.clone(), e.into()))?;
            let mut writer = FastqWriter::create(out1).map_err(|e| (lines.clone(), e.into()))?;
            while let Some(rec) = reader.next_record() {
                let rec = rec.map_err(|e| (lines.clone(), e.into()))?;
                seen += 1;
                if progress && seen.is_multiple_of(1_000_000) {
                    lines.info(format!("    Processed {seen} reads"));
                }
                if dups.contains_name(rec.name0()) {
                    removed += 1;
                } else {
                    writer
                        .write_record(&rec.name, &rec.seq, &rec.qual)
                        .map_err(|e| (lines.clone(), e.into()))?;
                    kept += 1;
                }
            }
            writer.finish().map_err(|e| (lines.clone(), e.into()))?;
        }
        // The output tuple mirrors the input tuple (M1), so a mixed
        // single/pair shape cannot occur.
        _ => unreachable!("output tuple mirrors input tuple"),
    }

    lines.info(format!("    Removed {removed} reads"));
    lines.info(format!("    Kept {kept} reads"));
    lines.info(format!(
        "    Purged fastq file saved to '{}'",
        out1.display()
    ));
    Ok(lines)
}

/// `file_name` that renders as the empty string rather than erroring on
/// paths without a file name (defensive; discovery always yields files).
pub(crate) fn file_name(p: &Path) -> String {
    p.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}
