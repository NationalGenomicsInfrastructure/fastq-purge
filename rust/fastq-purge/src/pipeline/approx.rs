//! The approx (bloom) pipeline (milestone M4), canonical semantics per
//! RUST_REWRITE.md §3/§5 (decisions D1–D3) plus decision D4 (§8/M4):
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

use crate::cli::AppConfig;
use crate::discovery::purged_name;
use crate::fastq::{FastqWriter, FileFastqReader};
use crate::membership::bloom::{hash_name0 as hash_name0_bloom, Bloom64};
use crate::membership::exact::{hash_name0, ExactSet};
use crate::pipeline::exact::{file_name, write_purged_undetermined};

pub fn run(config: &AppConfig) -> anyhow::Result<()> {
    config.log.info("Using bloom filter method");

    // ---- pass 1: one bloom per undetermined lane -------------------
    config.log.info("Building bloom filter...");
    let mut blooms: Vec<Bloom64> = Vec::new();
    for (_lane, (r1, _r2)) in &config.undetermined {
        config
            .log
            .debug(&format!("Adding '{}' to bloom filter", file_name(r1)));
        let mut bloom = Bloom64::new(config.max_items, config.fpr);
        let mut reader = FileFastqReader::open(r1)?;
        while let Some(name0) = reader.next_name0() {
            bloom.insert_name(name0?);
        }
        let n = bloom.observed();
        let mb = bloom.estimated_bytes() as f64 / 1024.0 / 1024.0;
        config
            .log
            .info(&format!("Number of undetermined reads: {n} ({mb:.2} MB)"));
        config.log.debug(&format!(
            "Bloom filter: {} bits, {} probes, {} items",
            bloom.bits(),
            bloom.k(),
            bloom.capacity()
        ));
        blooms.push(bloom);
    }

    if config.purge_assigned {
        return purge_assigned(config, &blooms);
    }

    // ---- pass 2: one scan of the assigned files, all lane blooms ----
    config
        .log
        .info("Scanning assigned files for already-assigned reads...");
    let k = blooms.len();
    let mut dups: Vec<ExactSet> = (0..k).map(|_| ExactSet::new()).collect();
    for (_lane, files) in &config.assigned {
        for f in files {
            config
                .log
                .info(&format!("Processing target file '{}'", file_name(f)));
            let mut reader = FileFastqReader::open(f)?;
            let mut new_here = 0usize;
            while let Some(name0) = reader.next_name0() {
                let name0 = name0?;
                let (h1, h2) = hash_name0_bloom(name0);
                let h64 = hash_name0(name0);
                for l in 0..k {
                    if blooms[l].contains_h(h1, h2) && dups[l].insert_hash(h64) {
                        new_here += 1;
                    }
                }
            }
            config
                .log
                .debug(&format!("'{}': {new_here} new duplicates", file_name(f)));
        }
    }

    // ---- pass 3: write the purged undetermined files -----------------
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
            true,
        )?;
    }

    config.log.info("Done!");
    Ok(())
}

/// D3: purge the ASSIGNED files instead (assigned reads whose name0
/// appears in the lane's undetermined bloom are removed). Assigned files
/// are matched to undetermined lanes by lane key; an assigned lane with
/// no undetermined lane is skipped with a warning. The membership test is
/// a bloom query, so a structural false positive here *can* remove a read
/// that should be kept (this is the only place the FPR is observable in
/// the outputs).
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
            let mut seen = 0u64;
            while let Some(rec) = reader.next_record() {
                let rec = rec?;
                seen += 1;
                if seen.is_multiple_of(1_000_000) {
                    config.log.info(&format!("    Processed {seen} reads"));
                }
                if blooms[l].contains_name(rec.name0()) {
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
