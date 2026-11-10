//! Dev/verification tool (M7.1): dump the parsed sample-sheet structure as
//! `lane<TAB>project<TAB>name` lines (one per sample, in the canonical
//! (Lane, Sample_Project, row) order). Used for differential testing
//! against the Python `parse_sample_sheet` (samplesheet_support branch).
//!
//! Usage: sheet_dump <sheet.csv>

fn main() {
    let path = std::env::args()
        .nth(1)
        .expect("usage: sheet_dump <sheet.csv>");
    match fastq_purge::samplesheet::parse(std::path::Path::new(&path)) {
        Ok(data) => {
            for (lane, projects) in &data.lanes {
                for (project, names) in projects {
                    for name in names {
                        println!("{lane}\t{project}\t{name}");
                    }
                }
            }
        }
        Err(e) => {
            eprintln!("{e}");
            std::process::exit(1);
        }
    }
}
