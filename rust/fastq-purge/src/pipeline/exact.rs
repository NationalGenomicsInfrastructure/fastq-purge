//! The exact pipeline (milestone M3), canonical semantics per
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
//! file lines; see `purge_assigned_lane`.)

use std::path::Path;

use crate::cli::AppConfig;
use crate::discovery::purged_name;
use crate::fastq::{FastqWriter, FileFastqReader, PairedFastqWriter};
use crate::membership::exact::{hash_name0, ExactSet};

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    config.log.info("Using exact method");

    // ---- pass 1: membership set per undetermined lane ----------------
    config
        .log
        .info("Building membership set from undetermined reads...");
    let mut sets: Vec<ExactSet> = Vec::new();
    for (_lane, (r1, _r2)) in &config.undetermined {
        config
            .log
            .debug(&format!("Adding '{}' to set", file_name(r1)));
        let mut set = ExactSet::new();
        let mut reader = FileFastqReader::open(r1)?;
        while let Some(name0) = reader.next_name0() {
            set.insert_name(name0?);
        }
        let n = set.len();
        let mb = set.estimated_bytes() as f64 / 1024.0 / 1024.0;
        config
            .log
            .info(&format!("Number of undetermined reads: {n} ({mb:.2} MB)"));
        sets.push(set);
    }

    if config.purge_assigned {
        return purge_assigned(config, &sets);
    }

    // ---- pass 2: one scan of the assigned files, all lane sets --------
    config
        .log
        .info("Scanning assigned files for already-assigned reads...");
    let k = sets.len();
    let mut dups: Vec<ExactSet> = (0..k).map(|_| ExactSet::new()).collect();
    for (_lane, files) in &config.assigned {
        for f in files {
            config
                .log
                .info(&format!("Processing target file '{}'", file_name(f)));
            let mut reader = FileFastqReader::open(f)?;
            let mut new_here = 0usize;
            while let Some(name0) = reader.next_name0() {
                let h = hash_name0(name0?);
                for l in 0..k {
                    if sets[l].contains_hash(h) && dups[l].insert_hash(h) {
                        new_here += 1;
                    }
                }
            }
            config
                .log
                .debug(&format!("'{}': {new_here} new duplicates", file_name(f)));
        }
    }

    // ---- pass 3: write the purged undetermined files -------------------
    for (l, (_lane, (r1, r2))) in config.undetermined.iter().enumerate() {
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
        let (_lane, (out1, out2)) = &config.output[l];
        write_purged_undetermined(
            r1,
            r2.as_deref(),
            out1,
            out2.as_deref(),
            &dups[l],
            &config.log,
            false,
        )?;
    }

    config.log.info("Done!");
    Ok(())
}

/// D3: purge the ASSIGNED files instead (assigned reads whose name0
/// appears in the lane's undetermined set are removed). Assigned files
/// are matched to undetermined lanes by lane key; an assigned lane with
/// no undetermined lane is skipped with a warning.
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

    for (l, (lane, _und)) in config.undetermined.iter().enumerate() {
        let Some((_key, files)) = config.assigned.iter().find(|(k, _)| k == lane) else {
            config
                .log
                .debug(&format!("lane {lane}: no assigned files to purge"));
            continue;
        };
        for f in files {
            let out = config.output_dir.join(purged_name(f));
            config
                .log
                .info(&format!("Processing target file '{}'", file_name(f)));
            config
                .log
                .info(&format!("    Output file: '{}'", out.display()));
            let mut reader = FileFastqReader::open(f)?;
            let mut writer = FastqWriter::create(&out)?;
            let (mut kept, mut removed) = (0u64, 0u64);
            while let Some(rec) = reader.next_record() {
                let rec = rec?;
                if sets[l].contains_name(rec.name0()) {
                    removed += 1;
                } else {
                    writer.write_record(&rec.name, &rec.seq, &rec.qual)?;
                    kept += 1;
                }
            }
            writer.finish()?;
            config.log.info(&format!("    Removed {removed} reads"));
            config.log.info(&format!("    Kept {kept} reads"));
            config.log.info(&format!(
                "    Purged fastq file saved to '{}'",
                out.display()
            ));
        }
    }

    config.log.info("Done!");
    Ok(())
}

/// Write the purged undetermined file(s) of one lane (the default D1
/// direction). Paired inputs are zipped: iteration stops at the shorter
/// of R1/R2 (Python `dnaio` paired iteration is a zip). `progress`
/// emits the per-1M `    Processed {n} reads` line (the approx method's
/// contract; the exact method does not print it — Python's exact path
/// doesn't either). Shared by the exact and approx pipelines (M4: the
/// dups set is exact in both — decision D4).
pub(crate) fn write_purged_undetermined(
    r1: &Path,
    r2: Option<&Path>,
    out1: &Path,
    out2: Option<&Path>,
    dups: &ExactSet,
    log: &crate::logging::Logger,
    progress: bool,
) -> anyhow::Result<()> {
    let (mut kept, mut removed) = (0u64, 0u64);
    let mut seen = 0u64;

    match (r2, out2) {
        (Some(r2), Some(out2)) => {
            log.info(&format!(
                "Writing purged fastq files '{}' and '{}'",
                file_name(out1),
                file_name(out2)
            ));
            let mut reader1 = FileFastqReader::open(r1)?;
            let mut reader2 = FileFastqReader::open(r2)?;
            let mut writer = PairedFastqWriter::create(out1, out2)?;
            while let Some(rec1) = reader1.next_record() {
                let rec1 = rec1?;
                // zip semantics: stop at the shorter of R1/R2
                let Some(rec2) = reader2.next_record() else {
                    break;
                };
                let rec2 = rec2?;
                seen += 1;
                if progress && seen.is_multiple_of(1_000_000) {
                    log.info(&format!("    Processed {seen} reads"));
                }
                if dups.contains_name(rec1.name0()) {
                    removed += 1;
                } else {
                    writer.write_pair(&rec1, &rec2)?;
                    kept += 1;
                }
            }
            writer.finish()?;
        }
        (None, None) => {
            log.info(&format!("Writing purged fastq file '{}'", file_name(out1)));
            let mut reader = FileFastqReader::open(r1)?;
            let mut writer = FastqWriter::create(out1)?;
            while let Some(rec) = reader.next_record() {
                let rec = rec?;
                seen += 1;
                if progress && seen.is_multiple_of(1_000_000) {
                    log.info(&format!("    Processed {seen} reads"));
                }
                if dups.contains_name(rec.name0()) {
                    removed += 1;
                } else {
                    writer.write_record(&rec.name, &rec.seq, &rec.qual)?;
                    kept += 1;
                }
            }
            writer.finish()?;
        }
        // The output tuple mirrors the input tuple (M1), so a mixed
        // single/pair shape cannot occur.
        _ => unreachable!("output tuple mirrors input tuple"),
    }

    log.info(&format!("    Removed {removed} reads"));
    log.info(&format!("    Kept {kept} reads"));
    log.info(&format!(
        "    Purged fastq file saved to '{}'",
        out1.display()
    ));
    Ok(())
}

/// `file_name` that renders as the empty string rather than erroring on
/// paths without a file name (defensive; discovery always yields files).
pub(crate) fn file_name(p: &Path) -> String {
    p.file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default()
}
