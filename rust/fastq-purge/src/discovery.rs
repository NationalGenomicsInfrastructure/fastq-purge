//! File discovery and validation rules, ported from
//! `fastq_purge.main.validate_args` (Python).
//!
//! Parity notes (behaviour preserved from the Python implementation,
//! quirks included; see RUST_REWRITE.md §7):
//!
//! * `explode_path`: directories are scanned for entries whose *name*
//!   contains `.fq` (pattern `*.fq*`) or `.fastq` (pattern `*.fastq*`),
//!   the two pattern results concatenated in that order. Matches include
//!   directories and symlinks (as `pathlib.Path.glob` does); recursive
//!   search does not descend into symlinked directories. A file path is
//!   accepted as-is, whatever its name.
//! * `create_paired_dict`: each undetermined file's name is normalised
//!   with `re.sub(r"_R[12]_", "_R[12]_", name)` — i.e. a `_R1_` or `_R2_`
//!   segment becomes the *literal* text `_R[12]_`. The (parent, pattern)
//!   pairs are deduplicated, sorted, and re-globbed: the `[12]` character
//!   class then re-matches the R1 and R2 files, so a lane with both reads
//!   yields the pair `(R1, R2)` (sorted) and a lane with one read yields
//!   `(file, None)`. The lane key is the 3rd `_`-separated token of the
//!   *pattern* string. Later patterns with the same key overwrite earlier
//!   ones (dict semantics), insertion order is preserved.
//! * `clean_assigned`: files are sorted and deduplicated; the collapse
//!   `re.sub(r"_[IR][0-9]_", "_", name)` (all occurrences) is the
//!   identity key — a second file collapsing to a seen key is skipped with
//!   the "paired or index file" warning.
//! * The lane key is the first `re.search(r"_L[0-9]{3}_")` match, stripped
//!   of underscores.
//!
//! Deliberate, documented deviations from the Python CLI:
//!
//! * The Python CLI raises a raw traceback (IndexError) when an
//!   undetermined directory yields no files after filtering; this code
//!   returns a clean error (exit 1).
//! * An undetermined file name with fewer than three `_`-separated tokens
//!   crashes the Python CLI (`pattern.split("_")[2]`); this code returns a
//!   clean error.
//! * Glob metacharacters in file names: fnmatch backslash-escaping is not
//!   reproduced (backslash is a literal); `^` is not a class negator
//!   (only `!`, as in `fnmatch`).

use std::collections::HashMap;
use std::path::{Path, PathBuf};

use crate::logging::Logger;

/// A lane's undetermined input: (R1 or the single file, R2 or None).
pub type UndeterminedLane = (PathBuf, Option<PathBuf>);

/// Port of `explode_path(path, recursive)`.
///
/// On a path that is neither file nor directory this logs the Python
/// error message and returns an error (the Python CLI then `exit(1)`s).
pub fn explode_path(path: &Path, recursive: bool, log: &Logger) -> Result<Vec<PathBuf>, String> {
    if path.is_dir() {
        let mut files = Vec::new();
        if recursive {
            collect_recursive(path, &mut files, |name| name.contains(".fq"));
            collect_recursive(path, &mut files, |name| name.contains(".fastq"));
        } else {
            collect_flat(path, &mut files, |name| name.contains(".fq"));
            collect_flat(path, &mut files, |name| name.contains(".fastq"));
        }
        Ok(files)
    } else if path.is_file() {
        Ok(vec![path.to_path_buf()])
    } else {
        let msg = format!(
            "Path '{}' is neither a file nor a directory",
            path.display()
        );
        log.error(&msg);
        Err(msg)
    }
}

fn entry_name(path: &Path) -> Option<&str> {
    path.file_name().and_then(|n| n.to_str())
}

fn collect_flat(dir: &Path, out: &mut Vec<PathBuf>, matches: fn(&str) -> bool) {
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let Some(name) = entry_name(&path) else {
            continue;
        };
        if matches(name) {
            out.push(path);
        }
    }
}

fn collect_recursive(dir: &Path, out: &mut Vec<PathBuf>, matches: fn(&str) -> bool) {
    collect_flat(dir, out, matches);
    let Ok(entries) = std::fs::read_dir(dir) else {
        return;
    };
    for entry in entries.flatten() {
        let p = entry.path();
        // Recurse into real subdirectories only (pathlib `**` does not
        // follow symlinked directories).
        if p.is_dir() && !p.is_symlink() {
            collect_recursive(&p, out, matches);
        }
    }
}

/// `re.sub(r"_R[12]_", "_R[12]_", name)`: replace every `_R1_` / `_R2_`
/// segment with the literal text `_R[12]_`.
fn normalize_read_segment(name: &str) -> String {
    let b = name.as_bytes();
    let mut out = String::with_capacity(name.len() + 2);
    let mut i = 0;
    while i < b.len() {
        if i + 3 < b.len()
            && b[i] == b'_'
            && b[i + 1] == b'R'
            && matches!(b[i + 2], b'1' | b'2')
            && b[i + 3] == b'_'
        {
            out.push_str("_R[12]_");
            i += 4;
        } else {
            out.push(b[i] as char);
            i += 1;
        }
    }
    out
}

/// Port of `create_paired_dict` (see module docs for the semantics).
///
/// Returns lanes in discovery (insertion) order.
pub fn create_paired_dict(files: Vec<PathBuf>) -> Result<Vec<(String, UndeterminedLane)>, String> {
    let files: Vec<PathBuf> = files
        .into_iter()
        .filter(|p| !entry_name(p).unwrap_or("").contains("purged"))
        .collect();

    // Deduplicated, sorted (parent, pattern) pairs.
    let mut patterns: Vec<(PathBuf, String)> = files
        .iter()
        .filter_map(|x| {
            entry_name(x).and_then(|name| {
                x.parent()
                    .map(|p| (p.to_path_buf(), normalize_read_segment(name)))
            })
        })
        .collect();
    patterns.sort();
    patterns.dedup();

    let mut index: HashMap<String, usize> = HashMap::new();
    let mut lanes: Vec<(String, UndeterminedLane)> = Vec::new();
    for (parent, pattern) in &patterns {
        // Parity with Python's `pattern.split("_")[2]`: names with fewer
        // than three underscore-separated tokens crash the Python CLI
        // (IndexError); we fail cleanly.
        let key = pattern
            .split('_')
            .nth(2)
            .ok_or_else(|| format!("Undetermined file pattern '{pattern}' has no lane token"))?
            .to_string();

        let matches = glob_component(parent, pattern)?;
        // Parity with Python: exactly two matches -> (first, second)
        // (sorted); any other count -> (first, None). Zero matches cannot
        // happen for a pattern derived from an existing file, unless the
        // original name carries other glob metacharacters; fail cleanly.
        let lane: UndeterminedLane = match matches.len() {
            2 => (matches[0].clone(), Some(matches[1].clone())),
            1 => (matches[0].clone(), None),
            _ => {
                return Err(format!(
                    "Pattern '{pattern}' in '{}' matched no files",
                    parent.display()
                ))
            }
        };

        match index.get(&key) {
            Some(&pos) => lanes[pos].1 = lane, // overwrite, keep position
            None => {
                index.insert(key.clone(), lanes.len());
                lanes.push((key, lane));
            }
        }
    }
    Ok(lanes)
}

/// Single-component glob under `dir` with fnmatch metacharacters
/// (`*`, `?`, `[...]`, `[!...]`), mirroring `pathlib.Path.glob` on a
/// pattern without `/`.
fn glob_component(dir: &Path, pattern: &str) -> Result<Vec<PathBuf>, String> {
    let mut out = Vec::new();
    let Ok(entries) = std::fs::read_dir(dir) else {
        return Ok(out);
    };
    for entry in entries.flatten() {
        let path = entry.path();
        let Some(name) = entry_name(&path) else {
            continue;
        };
        if fnmatch(name, pattern) {
            out.push(path);
        }
    }
    out.sort();
    Ok(out)
}

/// fnmatch-style matching of `name` against a single-component `pattern`
/// (`*`, `?`, `[...]`, `[!...]`; backslash is a literal).
pub fn fnmatch(name: &str, pattern: &str) -> bool {
    match_fn(name.as_bytes(), pattern.as_bytes())
}

fn match_fn(name: &[u8], pattern: &[u8]) -> bool {
    let (mut n, mut p) = (0usize, 0usize);
    // (pattern index of the last '*', name index to retry at)
    let mut star: Option<(usize, usize)> = None;

    while n < name.len() {
        // (matched, pattern bytes to advance past the match). A matched
        // char class advances past the whole `[...]` expression; a
        // `[` without a closing `]` is a literal (fnmatch).
        let (matched, advance) = if p < pattern.len() {
            match pattern[p] {
                b'*' => {
                    star = Some((p, n));
                    p += 1;
                    // Re-evaluate: a star matches empty and we must look
                    // at the next pattern byte against the same name byte.
                    continue;
                }
                b'[' => match char_class(pattern, p, name[n]) {
                    Some((len, ok)) => (ok, len),
                    None => (name[n] == b'[', 1),
                },
                b'?' => (true, 1),
                c => (c == name[n], 1),
            }
        } else {
            (false, 0)
        };
        if matched {
            n += 1;
            p += advance;
            continue;
        }
        // Mismatch: backtrack to the most recent star. The star must consume one
        // more name character on every retry (and the retry point is
        // remembered), or patterns like `*.fastq` against `a.b.fastq`
        // retry the same state forever.
        match star {
            Some((sp, sn)) => {
                let sn = sn + 1;
                star = Some((sp, sn));
                n = sn;
                p = sp + 1;
            }
            _ => return false,
        }
    }
    while p < pattern.len() && pattern[p] == b'*' {
        p += 1;
    }
    p == pattern.len()
}

/// `[...]` at `pattern[p]`: returns (expression length, whether `c` is in
/// the class). A `[` with no closing `]` is not a class (None).
fn char_class(pattern: &[u8], p: usize, c: u8) -> Option<(usize, bool)> {
    if p + 1 >= pattern.len() {
        return None;
    }
    let mut i = p + 1;
    let negated = pattern[i] == b'!';
    if negated {
        i += 1;
    }
    let members_start = i;
    // A ']' immediately after '[' (or '[!') is a literal member.
    let mut j = i;
    if j < pattern.len() && pattern[j] == b']' {
        j += 1;
    }
    while j < pattern.len() && pattern[j] != b']' {
        j += 1;
    }
    if j >= pattern.len() {
        return None; // no closing bracket
    }
    let end = j;
    let mut in_class = false;
    let mut k = members_start;
    while k < end {
        if k + 2 < end && pattern[k + 1] == b'-' {
            let (lo, hi) = (pattern[k], pattern[k + 2]);
            if lo <= hi {
                in_class |= (lo..=hi).contains(&c);
            }
            k += 3;
        } else {
            in_class |= c == pattern[k];
            k += 1;
        }
    }
    Some((end - p + 1, in_class ^ negated))
}

/// Output file name for a purged input (parity with the naming code in
/// `validate_args` and `process_target_file`):
/// `X.fastq.gz -> X.purged.fastq.gz`, `X.fastq -> X.purged.fastq`.
pub fn purged_name(input: &Path) -> String {
    let name = input
        .file_name()
        .unwrap_or_default()
        .to_string_lossy()
        .into_owned();
    if name.ends_with(".gz") {
        let stripped = &name[..name.len() - 3]; // drop ".gz"
        let fq_suffix = suffix(stripped); // e.g. ".fastq" ("" when none)
                                          // Python: input_file.with_suffix("").with_suffix(".purged{fq}.gz")
                                          // replaces the stripped name's last suffix, so the stem is kept.
        let stem = &stripped[..stripped.len() - fq_suffix.len()];
        format!("{stem}.purged{fq_suffix}.gz")
    } else {
        let fq_suffix = suffix(&name); // e.g. ".fastq"
        let stem = &name[..name.len() - fq_suffix.len()];
        format!("{stem}.purged{fq_suffix}")
    }
}

/// Python `Path.suffix` for a plain file name (the part after the last
/// dot, dot included; "" when the name has no dot or it is a dotfile).
fn suffix(name: &str) -> &str {
    match name.rfind('.') {
        Some(pos) if pos > 0 => &name[pos..],
        _ => "",
    }
}

/// `re.sub(r"_[IR][0-9]_", "_", name)`: collapse every `_I1_`-style
/// segment (single digit) to `_`, all occurrences, left to right.
pub fn collapse(name: &str) -> String {
    let b = name.as_bytes();
    let mut out = String::with_capacity(name.len());
    let mut i = 0;
    while i < b.len() {
        if i + 3 < b.len()
            && b[i] == b'_'
            && matches!(b[i + 1], b'I' | b'R')
            && b[i + 2].is_ascii_digit()
            && b[i + 3] == b'_'
        {
            out.push('_');
            i += 4;
        } else {
            out.push(b[i] as char);
            i += 1;
        }
    }
    out
}

/// First `re.search(r"_L[0-9]{3}_", name)` match, stripped of the
/// surrounding underscores (e.g. "L001").
pub fn lane_key(name: &str) -> Option<String> {
    let b = name.as_bytes();
    if b.len() < 6 {
        return None;
    }
    for i in 0..=b.len() - 6 {
        if b[i] == b'_'
            && b[i + 1] == b'L'
            && b[i + 2].is_ascii_digit()
            && b[i + 3].is_ascii_digit()
            && b[i + 4].is_ascii_digit()
            && b[i + 5] == b'_'
        {
            return Some(name[i + 1..i + 5].to_string());
        }
    }
    None
}

/// Port of the assigned-file cleaning loop (dedupe by collapsed name,
/// warn on skips). Returns (kept files in sorted order, count skipped).
pub fn clean_assigned(files: Vec<PathBuf>, log: &Logger) -> (Vec<PathBuf>, usize) {
    let mut seen: Vec<String> = Vec::new();
    let mut kept = Vec::with_capacity(files.len());
    let mut skipped = 0usize;
    for f in files {
        let Some(name) = entry_name(&f) else {
            continue;
        };
        let base = collapse(name);
        if seen.contains(&base) {
            log.warning(&format!(
                "Ignoring '{}' as it is either a paired or index file",
                f.display()
            ));
            skipped += 1;
        } else {
            seen.push(base);
            kept.push(f);
        }
    }
    (kept, skipped)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::fs;
    use std::sync::atomic::{AtomicUsize, Ordering};

    use crate::logging::Level;

    static TMP_COUNTER: AtomicUsize = AtomicUsize::new(0);

    struct TmpDir(PathBuf);

    impl TmpDir {
        fn new() -> Self {
            let n = TMP_COUNTER.fetch_add(1, Ordering::SeqCst);
            let p = std::env::temp_dir().join(format!("fqp_m1_{}_{}", std::process::id(), n));
            fs::create_dir_all(&p).unwrap();
            TmpDir(p)
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TmpDir {
        fn drop(&mut self) {
            let _ = fs::remove_dir_all(&self.0);
        }
    }

    fn touch(dir: &Path, name: &str) -> PathBuf {
        let p = dir.join(name);
        fs::write(&p, b"").unwrap();
        p
    }

    fn quiet() -> Logger {
        Logger::new(Level::Critical)
    }

    fn names(paths: &[PathBuf]) -> Vec<String> {
        paths
            .iter()
            .map(|p| p.file_name().unwrap().to_string_lossy().into_owned())
            .collect()
    }

    // ------------------------------------------------------------------
    // fnmatch
    // ------------------------------------------------------------------

    #[test]
    fn fnmatch_literals_and_stars() {
        assert!(fnmatch("a.fastq", "*.fastq"));
        assert!(fnmatch("a.b.fastq", "*.fastq"));
        assert!(!fnmatch("fastq", "*.fastq"));
        assert!(fnmatch("anything", "*"));
        assert!(fnmatch("abcfastq", "*fastq"));
        assert!(fnmatch("x_y", "*_*"));
        assert!(!fnmatch("x", "*_*"));
        assert!(fnmatch("a", "a"));
        assert!(!fnmatch("ab", "a"));
        assert!(!fnmatch("a", "ab"));
    }

    #[test]
    fn fnmatch_question_mark() {
        assert!(fnmatch("L001", "L??1"));
        assert!(!fnmatch("L01", "L??1"));
        assert!(fnmatch("a_R1_x", "a_?1_x"));
    }

    #[test]
    fn fnmatch_char_classes() {
        assert!(fnmatch("R1", "R[12]"));
        assert!(fnmatch("R2", "R[12]"));
        assert!(!fnmatch("R3", "R[12]"));
        // negation (fnmatch uses '!', not '^')
        assert!(fnmatch("R3", "R[!12]"));
        assert!(!fnmatch("R1", "R[!12]"));
        assert!(!fnmatch("R3", "R[^12]")); // '^' is a literal member
                                           // ranges
        assert!(fnmatch("R5", "R[1-9]"));
        assert!(!fnmatch("R10", "R[1-9]"));
        assert!(fnmatch("Ra", "R[a-z]"));
        assert!(!fnmatch("RA", "R[a-z]"));
        // ']' as first member is literal
        assert!(fnmatch("]", "[]a]"));
        assert!(fnmatch("a", "[]a]"));
        assert!(!fnmatch("b", "[]a]"));
        // unclosed '[' is a literal
        assert!(fnmatch("[a", "[a"));
        assert!(!fnmatch("[a", "x[a"));
        // the real-world pattern
        let pat = "Undetermined_S0_L001_R[12]_001.fastq.gz";
        assert!(fnmatch("Undetermined_S0_L001_R1_001.fastq.gz", pat));
        assert!(fnmatch("Undetermined_S0_L001_R2_001.fastq.gz", pat));
        assert!(!fnmatch("Undetermined_S0_L001_R3_001.fastq.gz", pat));
    }

    // ------------------------------------------------------------------
    // name transformations
    // ------------------------------------------------------------------

    #[test]
    fn normalize_read_segment_replaces_r1_r2() {
        assert_eq!(
            normalize_read_segment("Undetermined_S0_L001_R1_001.fastq.gz"),
            "Undetermined_S0_L001_R[12]_001.fastq.gz"
        );
        assert_eq!(
            normalize_read_segment("Undetermined_S0_L001_R2_001.fastq.gz"),
            "Undetermined_S0_L001_R[12]_001.fastq.gz"
        );
        assert_eq!(
            normalize_read_segment("S1_L001_R1_001.fastq.gz"),
            "S1_L001_R[12]_001.fastq.gz"
        );
        // multi-digit read numbers are untouched
        assert_eq!(normalize_read_segment("X_R10_Y"), "X_R10_Y");
        // re.sub consumes the match, so the shared underscore of adjacent
        // segments blocks the second one (Python parity)
        assert_eq!(normalize_read_segment("X_R1_R2_Y"), "X_R[12]_R2_Y");
        assert_eq!(normalize_read_segment("X_R1__R2_Y"), "X_R[12]__R[12]_Y");
        // no _R1_/_R2_ segment: unchanged
        assert_eq!(
            normalize_read_segment("Undetermined_S0_L001_001.fastq.gz"),
            "Undetermined_S0_L001_001.fastq.gz"
        );
    }

    #[test]
    fn collapse_removes_read_and_index_segments() {
        assert_eq!(collapse("S1_L001_R1_001.fastq.gz"), "S1_L001_001.fastq.gz");
        assert_eq!(collapse("S1_L001_R2_001.fastq.gz"), "S1_L001_001.fastq.gz");
        assert_eq!(collapse("X_I1_001.fastq.gz"), "X_001.fastq.gz");
        assert_eq!(collapse("X_I2_001.fastq.gz"), "X_001.fastq.gz");
        // multi-digit never collapses
        assert_eq!(
            collapse("S1_L001_R10_001.fastq.gz"),
            "S1_L001_R10_001.fastq.gz"
        );
        // re.sub consumes the matched segment including its trailing underscore,
        // so the shared underscore of adjacent segments blocks the second match
        // (Python parity): only the first of "_R1_R2_" is collapsed.
        assert_eq!(
            collapse("S1_L001_R1_R2_001.fastq.gz"),
            "S1_L001_R2_001.fastq.gz"
        );
        assert_eq!(
            collapse("S1_L001_I1_I2_001.fastq.gz"),
            "S1_L001_I2_001.fastq.gz"
        );
        // the two segments "_R1_" and "_R2_" of "_R1__R2_" each collapse to one
        // underscore (their boundary underscores are part of the matches)
        assert_eq!(
            collapse("S1_L001_R1__R2_001.fastq.gz"),
            "S1_L001__001.fastq.gz"
        );
        // untouched when absent
        assert_eq!(collapse("S1_L001_001.fastq.gz"), "S1_L001_001.fastq.gz");
    }

    #[test]
    fn lane_key_extracts_first_three_digit_lane() {
        assert_eq!(
            lane_key("Undetermined_S0_L001_R1_001.fastq.gz"),
            Some("L001".to_string())
        );
        assert_eq!(
            lane_key("S1_L042_R2_001.fastq.gz"),
            Some("L042".to_string())
        );
        assert_eq!(lane_key("no_lane.fastq"), None);
        assert_eq!(lane_key("L001"), None); // needs the surrounding underscores
        assert_eq!(lane_key("x_L12_y"), None); // needs exactly three digits
        assert_eq!(lane_key("a_L0123_b"), None); // four digits, no 3-digit match
                                                 // only matches with a leading underscore count
        assert_eq!(lane_key("L001_x_L002_y"), Some("L002".to_string()));
    }

    #[test]
    fn purged_name_inserts_purged_before_last_suffix() {
        assert_eq!(
            purged_name(Path::new("Undetermined_S0_L001_R1_001.fastq.gz")),
            "Undetermined_S0_L001_R1_001.purged.fastq.gz"
        );
        assert_eq!(purged_name(Path::new("X.fastq")), "X.purged.fastq");
        assert_eq!(purged_name(Path::new("X.fq.gz")), "X.purged.fq.gz");
        assert_eq!(purged_name(Path::new("X.fq")), "X.purged.fq");
        assert_eq!(purged_name(Path::new("noext")), "noext.purged");
        // dotfile: Python Path.suffix is "" -> whole name is the stem
        assert_eq!(purged_name(Path::new(".fastq")), ".fastq.purged");
        // with a parent directory, only the file name matters
        assert_eq!(
            purged_name(Path::new("/some/dir/X.fastq.gz")),
            "X.purged.fastq.gz"
        );
    }

    // ------------------------------------------------------------------
    // explode_path
    // ------------------------------------------------------------------

    #[test]
    fn explode_path_flat_keeps_pattern_concat_order_with_duplicates() {
        let d = TmpDir::new();
        // "a.fq.fastq.gz" matches BOTH patterns -> appears twice
        touch(d.path(), "a.fq.fastq.gz");
        touch(d.path(), "b.fastq");
        touch(d.path(), "c.fq.gz");
        touch(d.path(), "d.txt"); // no .fq/.fastq substring: never matched
        let files = explode_path(d.path(), false, &quiet()).unwrap();
        // Two pattern passes, concatenated, NO dedup (parity with the
        // Python glob concatenation). Within a pass the order is
        // directory order, so compare as a multiset. Note ".fastq" does
        // NOT contain the substring ".fq" — the patterns are disjoint
        // except for names like a.fq.fastq.gz.
        let mut got = names(&files);
        got.sort();
        assert_eq!(
            got,
            vec!["a.fq.fastq.gz", "a.fq.fastq.gz", "b.fastq", "c.fq.gz"]
        );
    }

    #[test]
    fn explode_path_recursive_descends_real_dirs_only() {
        let d = TmpDir::new();
        touch(d.path(), "top.fastq.gz");
        let sub = d.path().join("sub");
        fs::create_dir(&sub).unwrap();
        touch(&sub, "nested.fastq");
        touch(&sub, "nested.txt");
        let real = d.path().join("real");
        fs::create_dir(&real).unwrap();
        touch(&real, "deeper.fq.gz");
        // symlink to a directory must NOT be descended into (pathlib `**`)
        #[cfg(unix)]
        std::os::unix::fs::symlink(&sub, d.path().join("link")).unwrap();

        let files = explode_path(d.path(), true, &quiet()).unwrap();
        let names = names(&files);
        assert!(names.contains(&"top.fastq.gz".to_string()));
        assert!(names.contains(&"nested.fastq".to_string()));
        assert!(names.contains(&"deeper.fq.gz".to_string()));
        assert!(!names.iter().any(|n| n == "nested.txt" || n == "deeper.txt"));
        // the symlink itself matches by name? no: "link" has no .fq/.fastq
        assert!(!names.contains(&"link".to_string()));
    }

    #[test]
    fn explode_path_file_passthrough_ignores_name() {
        let d = TmpDir::new();
        let f = touch(d.path(), "weird_name.txt"); // no .fq/.fastq at all
        let files = explode_path(&f, false, &quiet()).unwrap();
        assert_eq!(files, vec![f]);
    }

    #[test]
    fn explode_path_neither_file_nor_dir_errors() {
        let d = TmpDir::new();
        let missing = d.path().join("does_not_exist");
        let err = explode_path(&missing, false, &quiet()).unwrap_err();
        assert_eq!(
            err,
            format!(
                "Path '{}' is neither a file nor a directory",
                missing.display()
            )
        );
    }

    // ------------------------------------------------------------------
    // create_paired_dict
    // ------------------------------------------------------------------

    #[test]
    fn paired_dict_pairs_r1_r2_and_handles_single_read() {
        let d = TmpDir::new();
        let r1 = touch(d.path(), "Undetermined_S0_L001_R1_001.fastq.gz");
        let r2 = touch(d.path(), "Undetermined_S0_L001_R2_001.fastq.gz");
        let only = touch(d.path(), "Undetermined_S0_L002_R1_001.fastq.gz");
        let files = explode_path(d.path(), false, &quiet()).unwrap();
        let lanes = create_paired_dict(files).unwrap();
        assert_eq!(lanes.len(), 2);
        assert_eq!(lanes[0].0, "L001");
        assert_eq!(lanes[0].1, (r1, Some(r2)));
        assert_eq!(lanes[1].0, "L002");
        assert_eq!(lanes[1].1, (only, None));
    }

    #[test]
    fn paired_dict_excludes_purged_files() {
        let d = TmpDir::new();
        let r1 = touch(d.path(), "Undetermined_S0_L001_R1_001.fastq.gz");
        let r2 = touch(d.path(), "Undetermined_S0_L001_R2_001.fastq.gz");
        // a previous run's output, still in the directory
        touch(d.path(), "Undetermined_S0_L001_R1_001.purged.fastq.gz");
        touch(d.path(), "Undetermined_S0_L001_R2_001.purged.fastq.gz");
        let files = explode_path(d.path(), false, &quiet()).unwrap();
        let lanes = create_paired_dict(files).unwrap();
        assert_eq!(lanes.len(), 1);
        assert_eq!(lanes[0].1, (r1, Some(r2)));
    }

    #[test]
    fn paired_dict_same_key_overwrite_keeps_first_position() {
        let d = TmpDir::new();
        // two runs (001, 002) of lane L001, plus lane L002
        for run in ["001", "002"] {
            touch(d.path(), &format!("Undetermined_S0_L001_R1_{run}.fastq.gz"));
            touch(d.path(), &format!("Undetermined_S0_L001_R2_{run}.fastq.gz"));
        }
        touch(d.path(), "Undetermined_S0_L002_R1_001.fastq.gz");
        touch(d.path(), "Undetermined_S0_L002_R2_001.fastq.gz");
        let files = explode_path(d.path(), false, &quiet()).unwrap();
        let lanes = create_paired_dict(files).unwrap();
        // L001 first (first pattern inserted), holding the LAST pattern's
        // glob results (dict overwrite); L002 second.
        assert_eq!(lanes.len(), 2);
        assert_eq!(lanes[0].0, "L001");
        assert_eq!(
            names(&[lanes[0].1 .0.clone(), lanes[0].1 .1.clone().unwrap()]),
            vec![
                "Undetermined_S0_L001_R1_002.fastq.gz",
                "Undetermined_S0_L001_R2_002.fastq.gz"
            ]
        );
        assert_eq!(lanes[1].0, "L002");
    }

    #[test]
    fn paired_dict_name_without_lane_token_errors() {
        let d = TmpDir::new();
        touch(d.path(), "x.fastq.gz"); // split("_")[2] would crash Python
        let files = explode_path(d.path(), false, &quiet()).unwrap();
        assert!(create_paired_dict(files).is_err());
    }

    // ------------------------------------------------------------------
    // clean_assigned
    // ------------------------------------------------------------------

    #[test]
    fn clean_assigned_dedupes_by_collapsed_name_with_warnings() {
        let d = TmpDir::new();
        touch(d.path(), "S1_L001_R1_001.fastq.gz");
        touch(d.path(), "S1_L001_R2_001.fastq.gz");
        let i1 = touch(d.path(), "S1_L001_I1_001.fastq.gz");
        touch(d.path(), "S1_L001_I2_001.fastq.gz");
        let r1_other = touch(d.path(), "S1_L002_R1_001.fastq.gz");
        let mut files: Vec<PathBuf> = explode_path(d.path(), true, &quiet()).unwrap();
        files.dedup();
        files.sort();
        let (kept, skipped) = clean_assigned(files, &quiet());
        // All four S1_L001 files collapse to "S1_L001_001.fastq.gz"; in
        // sorted order the I1 file comes first ("I" < "R"), so it is the
        // one that survives (Python does the same).
        assert_eq!(kept, vec![i1, r1_other]);
        assert_eq!(skipped, 3);
    }

    #[test]
    fn clean_assigned_keeps_distinct_collapsed_names() {
        let d = TmpDir::new();
        let a = touch(d.path(), "S1_L001_R1_001.fastq.gz");
        let b = touch(d.path(), "S1_L001_R1_002.fastq.gz");
        let mut files: Vec<PathBuf> = explode_path(d.path(), false, &quiet()).unwrap();
        files.dedup();
        files.sort();
        let (kept, skipped) = clean_assigned(files, &quiet());
        assert_eq!(kept, vec![a, b]);
        assert_eq!(skipped, 0);
    }
}
