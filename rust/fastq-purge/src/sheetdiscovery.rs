//! Sheet-mode file discovery: port of `process_args` from the
//! `samplesheet_support` branch (ss_main.py L319-552).
//!
//! Parity notes (behaviour preserved, quirks included):
//!
//! * Undetermined files: recursive explode, names must start with
//!   `Undetermined` (branch include-regex `^Undetermined`); `purged`
//!   exclusion happens inside `create_undetermined_paired_dict`.
//! * Assigned files: recursive explode, names containing `Undetermined`
//!   excluded; set-deduplicated, sorted, cleaned with the branch's
//!   `Ignoring '...' as another file with the same basename` DEBUG line
//!   (not master's WARNING message — different function in the branch).
//! * Lane keys come from the first `_L[0-9]{3}_` token (both undetermined
//!   and assigned).
//! * Selection iterates in sample-sheet lane/project order; assigned files
//!   match by sample-name prefix.
//!
//! Deliberate deviations (branch crashes, we fail/skip cleanly):
//!
//! * No `_L[0-9]{3}_` token in an undetermined pattern -> clean error
//!   (branch: AttributeError); same precedent as M1.
//! * No `SampleSheet.(csv|txt)` in the flowcell -> clean error (branch:
//!   IndexError on an empty list).
//! * Single-read `(file, None)` lanes in the `Lane '...' details:` block
//!   -> the `None` entry is skipped (branch: AttributeError on
//!   `None.name`, the f-string is evaluated eagerly before the skip gate).
//! * The branch's second no-lane check (`all(not assigned or not
//!   undetermined)`) is unreachable — every lane reaching it passed the
//!   per-lane gate — and is omitted.

use std::collections::BTreeSet;
use std::path::{Path, PathBuf};

use crate::discovery::{self, UndeterminedLane};
use crate::logging::Logger;
use crate::samplesheet::SheetData;

/// One sheet-mode lane, in sample-sheet order.
pub struct SheetLane {
    /// Lane key, e.g. "L001".
    pub key: String,
    /// (R1 or the single file, R2 or None) — a 2-element tuple in Python.
    pub undetermined: UndeterminedLane,
    /// (output R1, output R2 or None).
    pub output: (PathBuf, Option<PathBuf>),
    /// Assigned files, one inner list per sheet project (sheet order).
    /// Inner lists may be empty; the pipeline filters them out (parity
    /// with the branch's `[x for l in ... for x in l if x]`).
    pub assigned: Vec<Vec<PathBuf>>,
}

/// Validated sheet-mode plan.
pub struct SheetPlan {
    pub lanes: Vec<SheetLane>,
    pub output_dir: PathBuf,
}

fn fname(p: &Path) -> &str {
    p.file_name().and_then(|n| n.to_str()).unwrap_or("")
}

/// Port of the branch's `process_args` discovery/selection part.
///
/// Returns `Err((msg, already_logged))`: `already_logged == true` when the
/// message was already emitted by the logger (branch `logger.error` +
/// `exit(1)` sites).
pub fn plan(
    undetermined: &Path,
    assigned_paths: &[PathBuf],
    output: Option<&Path>,
    sheet: &SheetData,
    log: &Logger,
) -> Result<SheetPlan, (String, bool)> {
    // ---- undetermined ---------------------------------------------------
    let und_files = discovery::explode_path(undetermined, true, log).map_err(|e| (e, true))?;
    let und_files: Vec<PathBuf> = und_files
        .into_iter()
        .filter(|p| fname(p).starts_with("Undetermined"))
        .collect();

    let undetermined =
        discovery::create_undetermined_paired_dict(und_files).map_err(|e| (e, false))?;
    if undetermined.is_empty() {
        return Err((
            "No undetermined files found! Please provide valid paths.".to_string(),
            true,
        ));
    }

    log.debug("Undetermined files:");
    for (lane, (f1, f2)) in &undetermined {
        log.debug(&format!("    {lane}: '{}'", f1.display()));
        match f2 {
            Some(f) => log.debug(&format!("    {lane}: '{}'", f.display())),
            None => log.debug(&format!("    {lane}: 'None'")),
        }
    }

    // ---- output ----------------------------------------------------------
    let output_dir = match output {
        Some(p) => p.to_path_buf(),
        None => undetermined[0]
            .1
             .0
            .parent()
            .unwrap_or(Path::new("."))
            .to_path_buf(),
    };
    if !output_dir.is_dir() {
        log.debug("The output path does not exist! It will be created...");
        std::fs::create_dir_all(&output_dir).map_err(|e| {
            (
                format!(
                    "Failed to create output path '{}': {e}",
                    output_dir.display()
                ),
                false,
            )
        })?;
    }

    let output: Vec<(String, (PathBuf, Option<PathBuf>))> = undetermined
        .iter()
        .map(|(lane, (f1, f2))| {
            let o1 = output_dir.join(discovery::purged_name(f1));
            let o2 = f2
                .as_ref()
                .map(|f| output_dir.join(discovery::purged_name(f)));
            (lane.clone(), (o1, o2))
        })
        .collect();

    log.debug("Output files:");
    for (lane, (o1, o2)) in &output {
        log.debug(&format!("    {lane}: '{}'", o1.display()));
        match o2 {
            Some(o) => log.debug(&format!("    {lane}: '{}'", o.display())),
            None => log.debug(&format!("    {lane}: 'None'")),
        }
    }

    // ---- assigned --------------------------------------------------------
    let mut all: BTreeSet<PathBuf> = BTreeSet::new();
    for p in assigned_paths {
        let files = discovery::explode_path(p, true, log).map_err(|e| (e, true))?;
        for f in files {
            if fname(&f).contains("Undetermined") {
                continue;
            }
            all.insert(f);
        }
    }

    // Branch: sorted(set(...)); BTreeSet iteration is the same order.
    let mut seen: BTreeSet<String> = BTreeSet::new();
    let clean: Vec<PathBuf> = all
        .iter()
        .filter(|f| {
            let base = discovery::collapse(fname(f));
            if seen.insert(base.clone()) {
                true
            } else {
                log.debug(&format!(
                    "Ignoring '{}' as another file with the same basename ('{}') was found in the same path",
                    fname(f),
                    strip_fastq_suffix(&base)
                ));
                false
            }
        })
        .cloned()
        .collect();

    let mut assigned_map: Vec<(String, Vec<PathBuf>)> = Vec::new();
    for f in &clean {
        match discovery::lane_key(fname(f)) {
            Some(key) => match assigned_map.iter_mut().find(|(k, _)| *k == key) {
                Some((_, files)) => files.push(f.clone()),
                None => assigned_map.push((key, vec![f.clone()])),
            },
            None => log.warning(&format!(
                "File '{}' does not match the expected pattern '_L[0-9]{{3}}_'",
                f.display()
            )),
        }
    }

    log.debug("Assigned files:");
    for (key, files) in &assigned_map {
        for val in files {
            log.debug(&format!("    {key}: '{}'", val.display()));
        }
    }

    // ---- selection (sample-sheet lane/project order) ----------------------
    let mut lanes: Vec<SheetLane> = Vec::new();
    for (key, projects) in &sheet.lanes {
        let und = undetermined.iter().find(|(k, _)| k == key).map(|(_, v)| v);
        let out = output.iter().find(|(k, _)| k == key).map(|(_, v)| v);

        let assigned_lists: Vec<Vec<PathBuf>> = projects
            .iter()
            .map(|(_, names)| {
                assigned_map
                    .iter()
                    .find(|(k, _)| k == key)
                    .map(|(_, files)| {
                        files
                            .iter()
                            .filter(|f| names.iter().any(|n| fname(f).starts_with(n.as_str())))
                            .cloned()
                            .collect()
                    })
                    .unwrap_or_default()
            })
            .collect();

        if assigned_lists.is_empty() || und.is_none() {
            log.warning(&format!(
                "Lane '{key}' should undergo purging, but no undetermined or assigned reads were found! Skipping this lane."
            ));
            continue;
        }

        let (f1, f2) = und.unwrap();
        // The branch's tuple is always 2 elements ((f1, f2) or (f1, None)).
        log.info(&format!(
            "Lane '{key}' has 2 undetermined reads and {} assigned reads files",
            assigned_lists.len()
        ));

        // `None` entries are skipped (branch crashes on `None.name`);
        // trailing literal apostrophe is a branch typo, kept for parity.
        let mut purged_names = vec![fname(f1)];
        if let Some(f2f) = f2 {
            purged_names.push(fname(f2f));
        }
        let assigned_names: Vec<&str> = assigned_lists.iter().flatten().map(|f| fname(f)).collect();
        log.debug(&format!(
            "Lane '{key}' details:\n\
             The following files will be purged: {}\n\
             From reads in the following files: {}'",
            purged_names.join(", "),
            assigned_names.join(", ")
        ));

        lanes.push(SheetLane {
            key: key.clone(),
            undetermined: (f1.clone(), f2.clone()),
            // The output dict is built from the same keys as the
            // undetermined dict, so `out` is present whenever `und` is.
            output: out
                .cloned()
                .expect("output key present for every undetermined key"),
            assigned: assigned_lists,
        });
    }

    if lanes.is_empty() {
        return Err((
            "No lanes found for purging! Please provide valid paths or a sample sheet.".to_string(),
            true,
        ));
    }

    Ok(SheetPlan { lanes, output_dir })
}

/// Port of `re.sub(r"(_001)?.f(ast)?q(.gz)?", "", base)` — a single
/// leftmost match with Python's greedy optional groups (each group tries
/// the longer form first and backtracks to empty).
fn strip_fastq_suffix(base: &str) -> String {
    let b = base.as_bytes();
    for i in 0..b.len() {
        for &with_001 in &[true, false] {
            let mut j = i;
            if with_001 {
                if b[i..].starts_with(b"_001") {
                    j = i + 4;
                } else {
                    continue;
                }
            }
            // `.` (any char) then literal 'f'
            if j + 2 > b.len() || b[j + 1] != b'f' {
                continue;
            }
            let mut k = j + 2;
            for &with_ast in &[true, false] {
                if with_ast {
                    if k + 3 <= b.len() && &b[k..k + 3] == b"ast" {
                        k += 3;
                    } else {
                        continue;
                    }
                }
                if k < b.len() && b[k] == b'q' {
                    k += 1;
                    if b[k..].starts_with(b".gz") {
                        k += 3;
                    }
                    return format!("{}{}", &base[..i], &base[k..]);
                }
            }
        }
    }
    base.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn strip_fastq_suffix_matches_python_re() {
        assert_eq!(
            strip_fastq_suffix("M00150_S1_L001_R1_001.fastq"),
            "M00150_S1_L001_R1"
        );
        assert_eq!(strip_fastq_suffix("X_001.fastq.gz"), "X");
        assert_eq!(strip_fastq_suffix("X.fastq"), "X");
        assert_eq!(strip_fastq_suffix("X.fq.gz"), "X");
        assert_eq!(strip_fastq_suffix("nope"), "nope");
    }

    /// Small end-to-end plan() test: two undetermined reads for L001, a
    /// read pair that collapses to one basename, and a lane-less file.
    #[test]
    fn plan_groups_and_dedupes() {
        let dir = std::env::temp_dir().join(format!("fqp_sheetdisc_{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&dir);
        std::fs::create_dir_all(dir.join("run").join("Undetermined")).unwrap();
        std::fs::create_dir_all(dir.join("run").join("Assigned")).unwrap();
        for (sub, name) in [
            ("Undetermined", "Undetermined_S1_L001_R1_001.fastq.gz"),
            ("Undetermined", "Undetermined_S1_L001_R2_001.fastq.gz"),
            ("Assigned", "S1_L001_R1_001.fastq.gz"),
            ("Assigned", "S1_L001_R2_001.fastq.gz"),
            ("Assigned", "no_lane.fastq.gz"),
        ] {
            std::fs::write(dir.join("run").join(sub).join(name), "").unwrap();
        }

        let sheet = SheetData {
            lanes: vec![(
                "L001".to_string(),
                vec![("projA".to_string(), vec!["S1".to_string()])],
            )],
        };
        let log = Logger::new(crate::logging::Level::Debug);
        let plan = plan(&dir.join("run"), &[dir.join("run")], None, &sheet, &log).unwrap();
        assert_eq!(plan.lanes.len(), 1);
        let lane = &plan.lanes[0];
        assert_eq!(lane.key, "L001");
        // R1 kept, R2 collapsed onto the same basename; no_lane dropped
        // (lane-less WARNING).
        assert_eq!(lane.assigned.len(), 1);
        assert_eq!(lane.assigned[0].len(), 1);
        assert!(fname(&lane.assigned[0][0]).contains("R1"));
        assert_eq!(
            lane.output.0.file_name().unwrap().to_str().unwrap(),
            "Undetermined_S1_L001_R1_001.purged.fastq.gz"
        );
        let _ = std::fs::remove_dir_all(&dir);
    }
}
