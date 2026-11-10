//! The approx (bloom) pipeline (milestone M4; M5: rayon), canonical
//! semantics per RUST_REWRITE.md §3/§5 (decisions D1–D3) plus decision
//! D4 (§8/M4):
//!
//! - D1: the undetermined files are purged — an undetermined read is
//!   removed when its name0 appears in any assigned file.
//!   `--purge-assigned` (D3) inverts this.
//! - D2: the comparison key is name0.
//! - D4: pass 1 uses one `Bloom64` per undetermined lane (memory-bounded:
//!   §6 sizing formula, ~19.17 bits/name at `--fpr 1e-4` vs ~128 bits/name
//!   for the exact set); pass 2 tests each assigned name0 against all K
//!   lane blooms and collects the matches in an **exact** `ExactSet` of
//!   hashes (duplicates ≪ undetermined, so it stays small); pass 3
//!   reuses the exact pipeline's writer.
//!
//!   Consequence (proof in §8/M4): a name enters the dups set only if it
//!   actually appears in the assigned files, and an undetermined read is
//!   removed only if its name0 is in the dups set — so the purged
//!   outputs are **identical to the exact method for any `--fpr`**; the
//!   structural FPR of the pass-1 bloom (affecting only the `Found {K}`
//!   count, which may exceed the exact count by the number of false
//!   positives) can never change which reads are removed.
//!
//! M5 (rayon): as in the exact pipeline — workers return their log lines
//! and the orchestrator flushes them in canonical lane/file order, so
//! the log is byte-identical to the M4 sequential run for any
//! `--threads`.
//!
//! Log contract (INFO level; `--purge-assigned` variant noted below):
//! ```text
//! Using bloom filter method
//! Building bloom filter...
//!   (DEBUG) Adding '{und R1 file name}' to bloom filter     [per lane]
//!   (DEBUG) Bloom filter: {bits} bits, {k} probes, {n} items [per lane]
//!   Number of undetermined reads: {N} ({MB:.2} MB)          [per lane]
//! Scanning assigned files for already-assigned reads...
//!   Processing target file '{assigned file name}'           [per file]
//!   (DEBUG) '{assigned file name}': {k} new duplicates      [per file]
//! Found {K} duplicates in the undetermined file.            [per lane]
//! [lane with K > 0:]
//!   Writing purged fastq file '{out1 name}'
//!   | Writing purged fastq files '{out1 name}' and '{out2 name}'
//!   (every 1 M)    Processed {n} reads
//!       Removed {R} reads
//!       Kept {W} reads
//!       Purged fastq file saved to '{out path}'
//! [lane with K == 0:]
//!   Found no duplicates in the undetermined reads. Exiting...
//! Done!
//! ```
//! Deliberate deviations from the Python approx log: the Python path
//! built one shared bloom from all lanes and then *filtered the assigned
//! files* (the broken direction); here the blooms are per-lane and the
//! undetermined files are filtered (D1), so the Python
//! `Output file:` / `Reading target fastq file...` lines (per assigned
//! target) do not occur. `Number of undetermined reads` reports the
//! record count (the bloom does not retain an exact set) with the bloom
//! byte size; `Found {K}` may exceed the exact method's count by the
//! number of structural false positives (see the D4 proof above).

use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::cli::AppConfig;
use crate::discovery::purged_name;
use crate::fastq::{FastqWriter, FileFastqReader};
use crate::membership::bloom::{hash_name0 as hash_name0_bloom, Bloom64};
use crate::membership::exact::{hash_name0, ExactSet};
use crate::pipeline::exact::{file_name, write_purged_undetermined};
use crate::pipeline::{flush_lines, Block, Lines, ScanBlock};

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    config.log.info("Using bloom filter method");

    // ---- pass 1: one bloom per undetermined lane -------------------
    config.log.info("Building bloom filter...");
    let built: Vec<Result<(Bloom64, Lines), (Lines, anyhow::Error)>> = config
        .undetermined
        .par_iter()
        .map(|(_lane, (r1, _r2))| build_bloom(config.max_items, config.fpr, r1))
        .collect();
    let mut blooms: Vec<Bloom64> = Vec::with_capacity(built.len());
    for b in built {
        match b {
            Ok((bloom, lines)) => {
                flush_lines(&config.log, &lines.v);
                blooms.push(bloom);
            }
            Err((lines, e)) => {
                flush_lines(&config.log, &lines.v);
                return Err(e);
            }
        }
    }

    if config.purge_assigned {
        return purge_assigned(config, &blooms);
    }

    // ---- pass 2: one scan of the assigned files, all lane blooms ----
    config
        .log
        .info("Scanning assigned files for already-assigned reads...");
    let k = blooms.len();
    let files: Vec<&Path> = config
        .assigned
        .iter()
        .flat_map(|(_lane, fs)| fs.iter().map(PathBuf::as_path))
        .collect();
    let scanned: Vec<ScanBlock> = files.par_iter().map(|f| scan_file(f, &blooms)).collect();
    let mut dups: Vec<ExactSet> = (0..k).map(|_| ExactSet::new()).collect();
    for (i, s) in scanned.into_iter().enumerate() {
        let (mut parts, lines) = match s {
            Ok(v) => v,
            Err((lines, e)) => {
                flush_lines(&config.log, &lines.v);
                return Err(e);
            }
        };
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

    // ---- pass 3: write the purged undetermined files -----------------
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
                true,
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

/// Pass-1 worker: stream one undetermined R1 file into a `Bloom64`.
/// Returns the bloom plus the log lines in sequential emission order
/// (`Adding ... to bloom filter`, then the count line, then the
/// DEBUG size line).
pub(crate) fn build_bloom(
    max_items: u64,
    fpr: f64,
    r1: &Path,
) -> Result<(Bloom64, Lines), (Lines, anyhow::Error)> {
    let mut lines = Lines::default();
    lines.debug(format!("Adding '{}' to bloom filter", file_name(r1)));
    let mut bloom = Bloom64::new(max_items, fpr);
    let mut reader = FileFastqReader::open(r1).map_err(|e| (lines.clone(), e.into()))?;
    while let Some(name0) = reader.next_name0() {
        bloom.insert_name(name0.map_err(|e| (lines.clone(), e.into()))?);
    }
    let n = bloom.observed();
    let mb = bloom.estimated_bytes() as f64 / 1024.0 / 1024.0;
    lines.info(format!("Number of undetermined reads: {n} ({mb:.2} MB)"));
    lines.debug(format!(
        "Bloom filter: {} bits, {} probes, {} items",
        bloom.bits(),
        bloom.k(),
        bloom.capacity()
    ));
    Ok((bloom, lines))
}

/// Pass-2 worker: stream one assigned file; per name0, compute the
/// bloom probe hashes (once) and the exact hash, and collect the
/// per-lane partial dups sets (idempotent insertion — the new-duplicate
/// count is computed by the orchestrator at merge time).
pub(crate) fn scan_file(f: &Path, blooms: &[Bloom64]) -> ScanBlock {
    let k = blooms.len();
    let mut lines = Lines::default();
    lines.info(format!("Processing target file '{}'", file_name(f)));
    let mut parts: Vec<ExactSet> = (0..k).map(|_| ExactSet::new()).collect();
    let mut reader = FileFastqReader::open(f).map_err(|e| (lines.clone(), e.into()))?;
    while let Some(name0) = reader.next_name0() {
        let name0 = name0.map_err(|e| (lines.clone(), e.into()))?;
        let (h1, h2) = hash_name0_bloom(name0);
        let h64 = hash_name0(name0);
        for l in 0..k {
            if blooms[l].contains_h(h1, h2) {
                parts[l].insert_hash(h64);
            }
        }
    }
    Ok((parts, lines))
}

/// D3: purge the ASSIGNED files instead (assigned reads whose name0
/// appears in the lane's undetermined bloom are removed). Assigned files
/// are matched to undetermined lanes by lane key; an assigned lane with
/// no undetermined lane is skipped with a warning. The membership test
/// is a bloom query, so a structural false positive here *can* remove a
/// read that should be kept (this is the only place the FPR is
/// observable in the outputs). M5: per-file workers in parallel, blocks
/// flushed in lane/file order.
fn purge_assigned(config: &AppConfig, blooms: &[Bloom64]) -> anyhow::Result<()> {
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
                    .map(|f| purge_file(f, &config.output_dir, &blooms[l]))
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

/// `--purge-assigned` worker: one assigned file filtered against the
/// lane's undetermined bloom (progress line every 1 M reads, as in the
/// M4 sequential run).
fn purge_file(f: &Path, out_dir: &Path, bloom: &Bloom64) -> Block {
    let out = out_dir.join(purged_name(f));
    let mut lines = Lines::default();
    lines.info(format!("Processing target file '{}'", file_name(f)));
    lines.info(format!("    Output file: '{}'", out.display()));
    let mut reader = FileFastqReader::open(f).map_err(|e| (lines.clone(), e.into()))?;
    let mut writer = FastqWriter::create(&out).map_err(|e| (lines.clone(), e.into()))?;
    let (mut kept, mut removed) = (0u64, 0u64);
    let mut seen = 0u64;
    while let Some(rec) = reader.next_record() {
        let rec = rec.map_err(|e| (lines.clone(), e.into()))?;
        seen += 1;
        if seen.is_multiple_of(1_000_000) {
            lines.info(format!("    Processed {seen} reads"));
        }
        if bloom.contains_name(rec.name0()) {
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
