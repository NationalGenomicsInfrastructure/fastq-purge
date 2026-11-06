//! Entry point: parse -> validate -> dispatch.
//!
//! Milestone M2 (RUST_REWRITE.md §8): CLI + discovery/validation (M1) and
//! FASTQ streaming + name0 hashing (M2); the pipeline passes land in M3/M4.

use clap::Parser;

fn main() {
    let cli = fastq_purge::cli::Cli::parse();
    match fastq_purge::cli::validate(cli) {
        Ok(config) => match fastq_purge::pipeline::run(&config) {
            Ok(()) => {}
            Err(e) => {
                eprintln!("{e}");
                std::process::exit(1);
            }
        },
        Err((msg, code, already_logged)) => {
            if !already_logged {
                eprintln!("{msg}");
            }
            std::process::exit(code);
        }
    }
}
