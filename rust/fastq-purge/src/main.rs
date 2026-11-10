//! Entry point: parse -> validate -> dispatch.
//!
//! Milestone M2 (RUST_REWRITE.md §8): CLI + discovery/validation (M1) and
//! FASTQ streaming + name0 hashing (M2); the pipeline passes land in M3/M4.

use clap::Parser;

fn main() {
    let cli = fastq_purge::cli::Cli::parse();
    match fastq_purge::cli::validate(cli) {
        Ok(config) => {
            // M5: the pipeline passes run on the global rayon pool,
            // sized by --threads (RUST_REWRITE.md §5).
            if let Err(e) = rayon::ThreadPoolBuilder::new()
                .num_threads(config.threads.max(1) as usize)
                .build_global()
            {
                eprintln!("failed to build the thread pool: {e}");
                std::process::exit(1);
            }
            match fastq_purge::pipeline::run(&config) {
                Ok(()) => {}
                Err(e) => {
                    // An empty Display marks an error that was already
                    // reported through the logger (M7 sheet-mode write
                    // failure); printing it again would duplicate the line.
                    let msg = e.to_string();
                    if !msg.is_empty() {
                        eprintln!("{msg}");
                    }
                    std::process::exit(1);
                }
            }
        }
        Err((msg, code, already_logged)) => {
            if !already_logged {
                eprintln!("{msg}");
            }
            std::process::exit(code);
        }
    }
}
