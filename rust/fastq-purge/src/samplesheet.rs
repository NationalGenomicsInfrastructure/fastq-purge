//! SampleSheet.csv parsing (milestone M7.1), parity with the
//! `samplesheet_support` branch's `parse_sample_sheet()` (main.py L173–261,
//! see RUST_REWRITE.md §8/M7):
//!
//! - The CSV header row is the first line starting with `FCID` (Illumina
//!   sample sheets carry a preamble before it); a sheet without one is an
//!   error the CLI logs to stdout and exits 1 for.
//! - `Lane` → `zfill(3)` → `pad_start(4, 'L')` ("1" → "L001", "10" →
//!   "L010").
//! - `index`/`index2` contribute their character count (empty → 0) to the
//!   per-row `Recipe := Recipe-{len(index)}-{len(index2)}` string.
//! - Only lanes with ≥ 2 distinct concatenated Recipes survive (the branch
//!   filters with a polars semi-join on the per-lane Recipe n_unique).
//! - Rows are ordered by (Lane, Sample_Project) — a stable sort, so
//!   sample names within a (lane, project) keep sheet row order.
//!
//! The Python implementation parses the tail with `polars.read_csv`; this
//! port is a dependency-free CSV reader (RFC-4180-style quotes) that
//! produces the identical structure. It is differential-tested against the
//! Python parser (tests/samplesheet.rs). Deliberate deviations from the
//! Python traceback behaviour (clean errors, M1 precedent): non-UTF-8
//! content, a missing required column, and a structurally broken CSV.

use std::fs;
use std::path::Path;

/// A parsed sample sheet, in the branch's ordering contract:
/// lanes in (Lane, Sample_Project) order; per lane, projects in sort
/// order; per project, sample names in sheet row order.
/// Sample names of one project, in sheet-row order.
pub type ProjectSamples = Vec<String>;
/// (Sample_Project, its sample names) pairs, in (Lane, Sample_Project)
/// order within a lane.
pub type LaneProjects = Vec<(String, ProjectSamples)>;
/// (Lane, its projects) pairs, in (Lane, Sample_Project) order — the
/// lanes that have ≥2 distinct recipes.
pub type LaneSelection = Vec<(String, LaneProjects)>;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct SheetData {
    pub lanes: LaneSelection,
}

impl SheetData {
    /// Flattened sample names of one lane (in (project, row) order) — the
    /// assigned-file prefix filter's match set.
    pub fn names_flat(&self, lane: &str) -> Vec<String> {
        self.lanes
            .iter()
            .find(|(l, _)| l == lane)
            .map(|(_, projects)| {
                projects
                    .iter()
                    .flat_map(|(_, names)| names.iter().cloned())
                    .collect()
            })
            .unwrap_or_default()
    }
}

/// Parsing failure. `NoFcid` is the one error the Python CLI logs to
/// stdout and exits 1 for; everything else traceback'd in Python and
/// fails cleanly here.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SheetError {
    NoFcid,
    Other(String),
}

impl std::fmt::Display for SheetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SheetError::NoFcid => write!(
                f,
                "No valid SampleSheet found! 'FCID' header not found in the file. "
            ),
            SheetError::Other(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for SheetError {}

/// Parse a sample sheet file (the branch reads `SampleSheet.csv`; the
/// `.txt` variant is accepted upstream and is content-identical here).
pub fn parse(path: &Path) -> Result<SheetData, SheetError> {
    let text = fs::read_to_string(path)
        .map_err(|e| SheetError::Other(format!("cannot read '{}': {e}", path.display())))?;

    // The FCID anchor: first line starting with "FCID" (the preamble
    // above it is ignored).
    let lines: Vec<&str> = text.split('\n').collect();
    let Some(skip) = lines.iter().position(|l| l.starts_with("FCID")) else {
        return Err(SheetError::NoFcid);
    };
    let csv_lines: Vec<&str> = lines[skip..]
        .iter()
        .map(|l| l.strip_suffix('\r').unwrap_or(*l))
        .collect();

    let header = parse_csv_line(csv_lines[0])
        .map_err(|e| SheetError::Other(format!("cannot parse sample sheet header: {e}")))?;
    let idx = |col: &str| -> Result<usize, SheetError> {
        header
            .iter()
            .position(|h| h == col)
            .ok_or_else(|| SheetError::Other(format!("missing column '{col}' in sample sheet")))
    };
    let i_lane = idx("Lane")?;
    let i_recipe = idx("Recipe")?;
    let i_proj = idx("Sample_Project")?;
    let i_name = idx("Sample_Name")?;
    let i_idx1 = idx("index")?;
    let i_idx2 = idx("index2")?;

    // (lane, project, sample_name) in sheet row order, with the lane's
    // distinct concatenated Recipes collected as we go.
    let mut rows: Vec<(String, String, String)> = Vec::new();
    let mut lane_recipes: Vec<(String, Vec<String>)> = Vec::new(); // distinct, row order
    for (no, line) in csv_lines[1..].iter().enumerate() {
        if line.is_empty() {
            continue; // trailing newline / blank line
        }
        let fields = parse_csv_line(line).map_err(|e| {
            SheetError::Other(format!("cannot parse sample sheet line {}: {e}", no + 2))
        })?;
        let field = |i: usize| fields.get(i).cloned().unwrap_or_default();
        // Lane: zfill(3) then pad_start(4, 'L') — string operations, as in
        // the branch (a non-numeric Lane yields the same garbled key).
        let lane_raw = field(i_lane);
        let lane = pad_start(&zfill(&lane_raw, 3), 4, 'L');
        let i1 = field(i_idx1).chars().count();
        let i2 = field(i_idx2).chars().count();
        let recipe = format!("{}-{i1}-{i2}", field(i_recipe));
        match lane_recipes.iter_mut().find(|(l, _)| l == &lane) {
            Some((_, rs)) => {
                if !rs.contains(&recipe) {
                    rs.push(recipe);
                }
            }
            None => lane_recipes.push((lane.clone(), vec![recipe])),
        }
        rows.push((lane.clone(), field(i_proj), field(i_name)));
    }

    // Keep only lanes with >= 2 distinct concatenated Recipes (the
    // branch's per-lane n_unique > 1 semi-join).
    let kept: Vec<String> = lane_recipes
        .iter()
        .filter(|(_, rs)| rs.len() >= 2)
        .map(|(l, _)| l.clone())
        .collect();

    // Stable sort by (Lane, Sample_Project); names within a group keep
    // row order (parity with polars' stable .sort([Lane, Sample_Project])).
    let mut order: Vec<usize> = (0..rows.len()).collect();
    order.sort_by(|&a, &b| {
        rows[a]
            .0
            .cmp(&rows[b].0)
            .then_with(|| rows[a].1.cmp(&rows[b].1))
    });

    let mut data = SheetData::default();
    for &i in &order {
        let (lane, proj, name) = &rows[i];
        if !kept.contains(lane) {
            continue;
        }
        let pos = data.lanes.iter().position(|(l, _)| l == lane);
        let projects: &mut Vec<(String, Vec<String>)> = match pos {
            Some(p) => &mut data.lanes[p].1,
            None => {
                data.lanes.push((lane.clone(), Vec::new()));
                let last = data.lanes.len() - 1;
                &mut data.lanes[last].1
            }
        };
        match projects.iter_mut().find(|(p, _)| p == proj) {
            Some((_, names)) => names.push(name.clone()),
            None => projects.push((proj.clone(), vec![name.clone()])),
        }
    }
    Ok(data)
}

/// Split one CSV line into fields (RFC-4180-style: double-quoted fields
/// may contain commas/quotes; `""` escapes a quote).
fn parse_csv_line(line: &str) -> Result<Vec<String>, String> {
    let mut fields: Vec<String> = Vec::new();
    let mut cur = String::new();
    let mut in_quotes = false;
    let mut chars = line.chars().peekable();
    while let Some(c) = chars.next() {
        if in_quotes {
            if c == '"' {
                if chars.peek() == Some(&'"') {
                    chars.next();
                    cur.push('"');
                } else {
                    in_quotes = false;
                }
            } else {
                cur.push(c);
            }
        } else {
            match c {
                '"' => in_quotes = true,
                ',' => {
                    fields.push(std::mem::take(&mut cur));
                }
                _ => cur.push(c),
            }
        }
    }
    if in_quotes {
        return Err(format!("unterminated quoted field in '{line}'"));
    }
    fields.push(cur);
    Ok(fields)
}

/// `str.zfill(3)` for ASCII (the branch's Lane values are ASCII digits).
fn zfill(s: &str, width: usize) -> String {
    if s.len() < width {
        format!("{}{}", "0".repeat(width - s.len()), s)
    } else {
        s.to_string()
    }
}

/// `str.pad_start(4, 'L')` for ASCII.
fn pad_start(s: &str, width: usize, pad: char) -> String {
    if s.len() < width {
        format!("{}{}", pad.to_string().repeat(width - s.len()), s)
    } else {
        s.to_string()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use std::path::PathBuf;

    /// Write `content` to a per-process temp file and return its path.
    fn sheet_file(tag: &str, content: &str) -> PathBuf {
        let p = std::env::temp_dir().join(format!("fqp_m7_{}_{}.csv", tag, std::process::id()));
        let mut f = std::fs::File::create(&p).unwrap();
        f.write_all(content.as_bytes()).unwrap();
        p
    }

    /// Standard header: the FCID row doubles as the polars header.
    const HDR: &str = "FCID,Plate_ID,Sample_Project,Sample_Name,index,index2,Recipe,Lane";

    #[test]
    fn basic_two_recipes_kept_single_recipe_lane_dropped() {
        let p = sheet_file(
            "basic",
            &format!(
                "{HDR}\nM00150.1.001,1,projA,S1,ATCACG,,A,1\nM00150.1.001,1,projA,S2,TAGCTA,,B,1\nM00150.1.001,1,projB,S3,ATCACG,,A,2\n"
            ),
        );
        let data = parse(&p).unwrap();
        // Lane 1 has recipes A-6-0 and B-6-0 (kept); lane 2 has only A-6-0
        // (dropped — absent from the result, no warning in the parser).
        assert_eq!(
            data.lanes,
            vec![(
                "L001".to_string(),
                vec![(
                    "projA".to_string(),
                    vec!["S1".to_string(), "S2".to_string()]
                )]
            )]
        );
        assert_eq!(
            data.names_flat("L001"),
            vec!["S1".to_string(), "S2".to_string()]
        );
        assert!(data.names_flat("L002").is_empty());
    }

    #[test]
    fn fcid_anchor_not_first_line() {
        let p = sheet_file(
            "anchor",
            &format!(
                "Illumina[BaseStation]FlowcellDescription\nM00150.1.001\n\n{HDR}\nM00150.1.001,1,projA,S1,ATCACG,,A,1\nM00150.1.001,1,projA,S2,TAGCTA,,B,1\n"
            ),
        );
        let data = parse(&p).unwrap();
        assert_eq!(data.lanes.len(), 1);
        assert_eq!(data.lanes[0].0, "L001");
    }

    #[test]
    fn quoted_comma_field() {
        let p = sheet_file(
            "quoted",
            &format!(
                "{HDR}\nM00150.1.001,1,\"Smith, John\",S1,ATCACG,,A,1\nM00150.1.001,1,\"Smith, John\",S2,TAGCTA,,B,1\n"
            ),
        );
        let data = parse(&p).unwrap();
        assert_eq!(
            data.lanes,
            vec![(
                "L001".to_string(),
                vec![(
                    "Smith, John".to_string(),
                    vec!["S1".to_string(), "S2".to_string()]
                )]
            )]
        );
    }

    /// Concatenated-recipe parity: same raw Recipe with different index
    /// lengths are distinct recipes → the lane survives the ≥2 filter.
    #[test]
    fn recipe_concatenation_index_lengths() {
        let p = sheet_file(
            "recipes",
            &format!(
                "{HDR}\nM00150.1.001,1,projA,S1,AAAA,,A,1\nM00150.1.001,1,projA,S2,AAAAAA,,A,1\n"
            ),
        );
        let data = parse(&p).unwrap();
        // A-4-0 vs A-6-0: two distinct recipes → kept.
        assert_eq!(data.lanes.len(), 1);
        // Empty index2 counts as 0; same index lengths would collapse.
        let p2 = sheet_file(
            "recipes_same",
            &format!(
                "{HDR}\nM00150.1.001,1,projA,S1,AAAA,,A,1\nM00150.1.001,1,projA,S2,TTTT,,A,1\n"
            ),
        );
        let data2 = parse(&p2).unwrap();
        assert!(
            data2.lanes.is_empty(),
            "A-4-0 twice must collapse to one recipe"
        );
    }

    #[test]
    fn lane_zfill_and_pad() {
        let p = sheet_file(
            "pad",
            &format!(
                "{HDR}\nM00150.1.001,1,p,S1,AAAA,,A,1\nM00150.1.001,1,p,S2,TTTT,,B,1\nM00150.1.001,1,p,S3,CCCC,,A,10\nM00150.1.001,1,p,S4,GGGG,,B,10\nM00150.1.001,1,p,S5,AAAA,,A,7\nM00150.1.001,1,p,S6,TTTT,,B,7\n"
            ),
        );
        let data = parse(&p).unwrap();
        let lanes: Vec<&str> = data.lanes.iter().map(|(l, _)| l.as_str()).collect();
        assert_eq!(lanes, vec!["L001", "L007", "L010"]);
    }

    #[test]
    fn stable_order_lane_project_row_order() {
        // Rows deliberately out of (lane, project) order; names within a
        // (lane, project) group keep sheet row order (stable sort).
        let p = sheet_file(
            "order",
            &format!(
                "{HDR}\nM00150.1.001,1,projB,P2b,AAAAAA,,A,1\nM00150.1.001,1,projA,P1a,AAAAAA,,B,1\nM00150.1.001,1,projA,P1a2,AAAAAA,,A,1\nM00150.1.001,1,projB,P2a,AAAAAA,,B,2\nM00150.1.001,1,projA,P1b,AAAAAA,,A,2\n"
            ),
        );
        let data = parse(&p).unwrap();
        assert_eq!(
            data.lanes,
            vec![
                (
                    "L001".to_string(),
                    vec![
                        (
                            "projA".to_string(),
                            vec!["P1a".to_string(), "P1a2".to_string()]
                        ),
                        ("projB".to_string(), vec!["P2b".to_string()]),
                    ]
                ),
                (
                    "L002".to_string(),
                    vec![
                        ("projA".to_string(), vec!["P1b".to_string()]),
                        ("projB".to_string(), vec!["P2a".to_string()]),
                    ]
                ),
            ]
        );
    }

    #[test]
    fn no_fcid_error_and_message() {
        let p = sheet_file("nofcid", "FCID,Plate_ID\nM00150,1\n");
        // Has an FCID header — but no required columns (exercises the
        // missing-column path; see no_fcid_below for the anchor case).
        let e = parse(&p).unwrap_err();
        assert_eq!(
            e,
            SheetError::Other("missing column 'Lane' in sample sheet".to_string())
        );

        let p2 = sheet_file("noanchor", "Illumina[BaseStation]\nM00150\n");
        let e2 = parse(&p2).unwrap_err();
        assert_eq!(e2, SheetError::NoFcid);
        // Exact branch message, trailing space included.
        assert_eq!(
            format!("{e2}"),
            "No valid SampleSheet found! 'FCID' header not found in the file. "
        );
    }

    #[test]
    fn missing_column_reports_first_in_check_order() {
        let p = sheet_file(
            "missing",
            "FCID,Plate_ID,Sample_Project,Sample_Name,index,index2,Lane\nM00150.1.001,1,p,S1,AAAA,,1\n",
        );
        let e = parse(&p).unwrap_err();
        assert_eq!(
            e,
            SheetError::Other("missing column 'Recipe' in sample sheet".to_string())
        );
    }

    #[test]
    fn crlf_and_blank_lines() {
        let p = sheet_file(
            "crlf",
            &format!(
                "{HDR}\r\nM00150.1.001,1,projA,S1,ATCACG,,A,1\r\n\r\nM00150.1.001,1,projA,S2,TAGCTA,,B,1\r\n"
            ),
        );
        let data = parse(&p).unwrap();
        assert_eq!(
            data.lanes,
            vec![(
                "L001".to_string(),
                vec![(
                    "projA".to_string(),
                    vec!["S1".to_string(), "S2".to_string()]
                )]
            )]
        );
    }
}
