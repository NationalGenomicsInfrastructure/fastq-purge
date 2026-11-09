//! Dev/verification tool for differential testing against the Python
//! implementation (M2, and M6's differential+perf validation):
//!
//! ```text
//! fastq-dump names   <in.fastq[.gz]> <out.txt>        one name0 per line
//! fastq-dump records <in.fastq[.gz]> <out.fastq[.gz]> normalised copy
//! ```
//!
//! `names` mirrors `record.name.split(" ")[0]` for every record; `records`
//! writes every record in the normalised `@name\nseq\n+\nqual\n` form
//! (gzip when the output extension is `.gz`).

use std::io::Write;
use std::path::Path;
use std::process::ExitCode;

use fastq_purge::fastq::{FastqWriter, FileFastqReader};

fn main() -> ExitCode {
    let args: Vec<String> = std::env::args().collect();
    if args.len() != 4 {
        eprintln!("usage: fastq-dump <names|records> <input.fastq[.gz]> <output.fastq[.gz]|.txt>");
        return ExitCode::from(2);
    }
    let result = match args[1].as_str() {
        "names" => dump_names(Path::new(&args[2]), Path::new(&args[3])),
        "records" => copy_records(Path::new(&args[2]), Path::new(&args[3])),
        other => {
            eprintln!("unknown mode '{other}' (expected 'names' or 'records')");
            return ExitCode::from(2);
        }
    };
    match result {
        Ok(n) => {
            println!("{n}");
            ExitCode::SUCCESS
        }
        Err(e) => {
            eprintln!("{e}");
            ExitCode::from(1)
        }
    }
}

fn dump_names(input: &Path, output: &Path) -> Result<u64, Box<dyn std::error::Error>> {
    let mut reader = FileFastqReader::open(input)?;
    let file = std::fs::File::create(output)?;
    let mut out = std::io::BufWriter::new(file);
    let mut n: u64 = 0;
    while let Some(rec) = reader.next_name0() {
        let name0 = rec?;
        out.write_all(name0)?;
        out.write_all(b"\n")?;
        n += 1;
    }
    out.flush()?;
    Ok(n)
}

fn copy_records(input: &Path, output: &Path) -> Result<u64, Box<dyn std::error::Error>> {
    let mut reader = FileFastqReader::open(input)?;
    let mut writer = FastqWriter::create(output)?;
    let mut n: u64 = 0;
    while let Some(rec) = reader.next_record() {
        let rec = rec?;
        writer.write_record_bytes(&rec.normalized_bytes())?;
        n += 1;
    }
    writer.finish()?;
    Ok(n)
}
