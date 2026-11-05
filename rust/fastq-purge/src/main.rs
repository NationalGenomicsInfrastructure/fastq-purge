//! Entry point: parse -> validate -> dispatch.
//!
//! Milestone M1 (RUST_REWRITE.md §8): CLI + discovery/validation port;
//! the pipeline passes land in M3/M4.

mod cli;
pub mod discovery;
mod logging;
mod pipeline;

use clap::Parser;

fn main() {
    let cli = cli::Cli::parse();
    match cli::validate(cli) {
        Ok(config) => match pipeline::run(&config) {
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
