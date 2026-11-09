//! FASTQ output writers.
//!
//! The Python side writes through dnaio with `compression_level=7` for
//! `.gz` outputs; records are always written in the normalised form
//! `@name\nseq\n+\nqual\n` (bare `+` line, LF terminators) regardless of
//! the input's line endings or `+`-line content.
//!
//! Note on `.gz` byte parity: Python's gzip stream (xopen) writes mtime=0,
//! but zlib and miniz-oxide at the same level do not produce identical
//! compressed bytes, so `.gz` outputs are compared after decompression.
//! Plain `.fastq` outputs are byte-identical.

use std::fs::File;
use std::io::{self, Write};
use std::path::Path;

use flate2::write::GzEncoder;
use flate2::Compression;

use super::read::Record;

/// Gzip level used by the Python implementation.
pub const GZIP_LEVEL: u32 = 7;

enum Inner {
    Plain(io::BufWriter<File>),
    Gz(GzEncoder<io::BufWriter<File>>),
}

impl Write for Inner {
    fn write(&mut self, buf: &[u8]) -> io::Result<usize> {
        match self {
            Inner::Plain(w) => w.write(buf),
            Inner::Gz(w) => w.write(buf),
        }
    }

    fn flush(&mut self) -> io::Result<()> {
        match self {
            Inner::Plain(w) => w.flush(),
            Inner::Gz(w) => w.flush(),
        }
    }
}

/// A single FASTQ output stream (plain or gzip, decided by extension).
pub struct FastqWriter {
    inner: Inner,
}

impl FastqWriter {
    pub fn create(path: &Path) -> io::Result<Self> {
        let file = File::create(path)?;
        let is_gz = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e == "gz");
        if is_gz {
            let enc = GzEncoder::new(io::BufWriter::new(file), Compression::new(GZIP_LEVEL));
            Ok(Self {
                inner: Inner::Gz(enc),
            })
        } else {
            Ok(Self {
                inner: Inner::Plain(io::BufWriter::new(file)),
            })
        }
    }

    /// Write one record as `@name\nseq\n+\nqual\n`.
    pub fn write_record(&mut self, name: &[u8], seq: &[u8], qual: &[u8]) -> io::Result<()> {
        self.inner.write_all(b"@")?;
        self.inner.write_all(name)?;
        self.inner.write_all(b"\n")?;
        self.inner.write_all(seq)?;
        self.inner.write_all(b"\n+\n")?;
        self.inner.write_all(qual)?;
        self.inner.write_all(b"\n")?;
        Ok(())
    }

    /// Write a pre-normalised record (`Record::normalized_bytes`).
    pub fn write_record_bytes(&mut self, rec: &[u8]) -> io::Result<()> {
        self.inner.write_all(rec)
    }

    /// Flush and finalise (gzip trailer), then sync to disk.
    pub fn finish(self) -> io::Result<()> {
        match self.inner {
            Inner::Plain(bw) => {
                bw.into_inner()?.sync_data()?;
            }
            Inner::Gz(enc) => {
                let bw = enc.finish()?;
                bw.into_inner()?.sync_data()?;
            }
        }
        Ok(())
    }
}

/// Two output streams written in lockstep (paired R1/R2 purging).
pub struct PairedFastqWriter {
    pub r1: FastqWriter,
    pub r2: FastqWriter,
}

impl PairedFastqWriter {
    pub fn create(path1: &Path, path2: &Path) -> io::Result<Self> {
        Ok(Self {
            r1: FastqWriter::create(path1)?,
            r2: FastqWriter::create(path2)?,
        })
    }

    pub fn write_pair(&mut self, r1: &Record, r2: &Record) -> io::Result<()> {
        self.r1.write_record(&r1.name, &r1.seq, &r1.qual)?;
        self.r2.write_record(&r2.name, &r2.seq, &r2.qual)?;
        Ok(())
    }

    pub fn finish(self) -> io::Result<()> {
        self.r1.finish()?;
        self.r2.finish()?;
        Ok(())
    }
}
