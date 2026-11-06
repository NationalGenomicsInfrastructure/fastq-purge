//! Streaming FASTQ(.gz) reader.
//!
//! Line-based port of dnaio 1.2.3's C parser; every error message and line
//! number below was verified empirically against it (RUST_REWRITE.md §8,
//! M2 findings). The reader works on a growing byte buffer: `next_name0`
//! borrows from the buffer (no per-record allocation on the membership
//! path), and `next_record` returns owned bytes for the write pass.

use std::fmt::{self, Display};
use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use flate2::read::MultiGzDecoder;

/// Chunk size for reads from the (decompressed) input stream.
const CHUNK: usize = 256 * 1024;
/// Re-pack the buffer once this many consumed bytes have accumulated.
const REPACK_AT: usize = 1024 * 1024;

/// A FASTQ format or I/O error. `Format` carries the exact
/// dnaio-compatible message.
#[derive(Debug, Clone)]
pub enum FastqError {
    Io(String),
    Format(String),
}

impl Display for FastqError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            FastqError::Io(msg) => f.write_str(msg),
            FastqError::Format(msg) => f.write_str(msg),
        }
    }
}

impl std::error::Error for FastqError {}

impl From<io::Error> for FastqError {
    fn from(e: io::Error) -> Self {
        FastqError::Io(e.to_string())
    }
}

fn format_error(line: Option<u64>, msg: String) -> FastqError {
    let where_ = match line {
        Some(n) => format!("line {n}"),
        None => "unknown line".to_string(),
    };
    FastqError::Format(format!("Error in FASTQ file at {where_}: {msg}"))
}

/// Python `repr()` of an ASCII byte string: single-quoted unless the string
/// contains a single quote (and no double quote), with backslash escapes for
/// the common controls. Used for the incomplete-record quote in
/// premature-EOF errors (the ASCII check guarantees the bytes are ASCII).
fn py_quote(s: &[u8]) -> String {
    let use_dq = s.contains(&b'\'') && !s.contains(&b'"');
    let q = if use_dq { '"' } else { '\'' };
    let mut out = String::with_capacity(s.len() + 2);
    out.push(q);
    for &b in s {
        match b {
            b'\n' => out.push_str("\\n"),
            b'\r' => out.push_str("\\r"),
            b'\t' => out.push_str("\\t"),
            b'\\' => out.push_str("\\\\"),
            b'\'' if !use_dq => out.push_str("\\'"),
            b'"' if use_dq => out.push_str("\\\""),
            c if c < 0x20 || c == 0x7f => out.push_str(&format!("\\x{c:02x}")),
            c => out.push(c as char),
        }
    }
    out.push(q);
    out
}

/// One parsed record. `name` is the header without the leading `@`/`>`
/// (line terminators normalised to LF, i.e. no trailing `\r`).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Record {
    pub name: Vec<u8>,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
}

impl Record {
    /// name0: the name up to the first space byte (0x20). Tabs are not
    /// separators; if there is no space, name0 is the whole name
    /// (`record.name.split(" ")[0]` in the Python code).
    pub fn name0(&self) -> &[u8] {
        match self.name.iter().position(|b| *b == b' ') {
            Some(i) => &self.name[..i],
            None => &self.name,
        }
    }

    /// The normalised 4-line form the writer emits: `@name\nseq\n+\nqual\n`
    /// (the `+` line is always bare, CRLF already normalised on read).
    pub fn normalized_bytes(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(self.name.len() + self.seq.len() + self.qual.len() + 4);
        out.push(b'@');
        out.extend_from_slice(&self.name);
        out.push(b'\n');
        out.extend_from_slice(&self.seq);
        out.push(b'\n');
        out.push(b'+');
        out.push(b'\n');
        out.extend_from_slice(&self.qual);
        out.push(b'\n');
        out
    }
}

/// Offsets of one parsed record inside the reader buffer, plus the name0
/// span. Line numbers are 1-based physical line numbers, as reported in
/// error messages.
#[derive(Debug)]
struct Line {
    no: u64,
    start: usize,
    end: usize,
}

/// A FASTQ stream reader over any `Read` (raw bytes or decompressed).
pub struct FastqReader<R: Read> {
    inner: R,
    /// Decompressed bytes; valid data is `buf[start..end)`.
    buf: Vec<u8>,
    start: usize,
    end: usize,
    /// Number of physical lines fully consumed so far.
    line_no: u64,
    /// Raw bytes of the current record attempt (with terminators as they
    /// appear in the stream), for the incomplete-record quote.
    rec: Vec<u8>,
}

impl<R: Read> FastqReader<R> {
    pub fn new(inner: R) -> Self {
        Self {
            inner,
            buf: Vec::new(),
            start: 0,
            end: 0,
            line_no: 0,
            rec: Vec::new(),
        }
    }

    /// Read one physical line. Returns its (1-based line number, start,
    /// end) in `self.buf`, with the `\n` and any trailing `\r` removed.
    /// A final unterminated line is returned once, then `None`.
    ///
    /// `allow_repack`: the line spans returned by earlier calls remain
    /// valid only while the buffer is not shifted, so repacking (which
    /// memmoves the buffered data to the front) is only permitted between
    /// records, never while a record's earlier line spans are still in
    /// use. Growth via `reserve` is fine at any point (no shift).
    fn read_line(&mut self, allow_repack: bool) -> Result<Option<Line>, FastqError> {
        loop {
            if let Some(rel) = self.buf[self.start..self.end]
                .iter()
                .position(|b| *b == b'\n')
            {
                let nl = self.start + rel;
                let mut line = Line {
                    no: self.line_no + 1,
                    start: self.start,
                    end: nl,
                };
                let mut content_end = nl;
                if content_end > line.start && self.buf[content_end - 1] == b'\r' {
                    content_end -= 1;
                }
                line.end = content_end;
                self.rec.extend_from_slice(&self.buf[line.start..nl + 1]);
                self.start = nl + 1;
                self.line_no = line.no;
                return Ok(Some(line));
            }
            // No newline in the buffered data: refill.
            if allow_repack && self.start >= REPACK_AT {
                let rem = self.end - self.start;
                self.buf.copy_within(self.start..self.end, 0);
                self.start = 0;
                self.end = rem;
            }
            if self.end == self.buf.len() {
                self.buf.reserve(CHUNK);
            }
            let pos = self.end;
            self.buf.resize(self.buf.capacity(), 0);
            let n = self.inner.read(&mut self.buf[pos..])?;
            self.end = pos + n;
            if n == 0 {
                if self.start < self.end {
                    // final line without a terminator
                    let mut line = Line {
                        no: self.line_no + 1,
                        start: self.start,
                        end: self.end,
                    };
                    if line.end > line.start && self.buf[line.end - 1] == b'\r' {
                        line.end -= 1;
                    }
                    self.rec.extend_from_slice(&self.buf[line.start..]);
                    self.start = self.end;
                    self.line_no = line.no;
                    return Ok(Some(line));
                }
                return Ok(None);
            }
        }
    }

    fn line_bytes(&self, line: &Line) -> &[u8] {
        &self.buf[line.start..line.end]
    }

    /// The `found '{c}'` character for a structural error on `line`: the
    /// first byte of the line, or the two characters `\n` for an empty line
    /// (as Python's f-string renders it).
    fn found_char(line: &Line, buf: &[u8]) -> String {
        if line.start == line.end {
            "\\n".to_string()
        } else {
            (buf[line.start] as char).to_string()
        }
    }

    /// Check one line for non-ASCII bytes (dnaio aborts as soon as the
    /// offending line is read, reported at "unknown line").
    fn ascii_check(&self, line: &Line) -> Result<(), FastqError> {
        if self.line_bytes(line).iter().any(|b| *b >= 0x80) {
            return Err(format_error(
                None,
                "Non-ASCII characters found in record.".to_string(),
            ));
        }
        Ok(())
    }

    /// The premature-EOF error: the next physical line number plus the raw
    /// bytes of the incomplete record attempt so far, Python-repr quoted.
    fn premature_eof(&self) -> FastqError {
        format_error(
            Some(self.line_no + 1),
            format!(
                "Premature end of file encountered. The incomplete final record was: {}",
                py_quote(&self.rec)
            ),
        )
    }

    /// Parse one full record. Returns the line spans; `None` at clean EOF
    /// before a new record started.
    ///
    /// Check order matches the C parser (empirically verified):
    /// 1. header read -> ASCII check (the `@` check is DEFERRED)
    /// 2. seq read (EOF -> premature) -> ASCII check
    /// 3. header must start with `@` or `>` (error at the header line)
    /// 4. plus line read (EOF -> premature) -> ASCII check
    /// 5. plus line must start with `+` (error at the plus line)
    /// 6. qual read (EOF -> premature) -> ASCII check
    /// 7. seq/qual lengths must match (error at the qual line)
    fn parse_record(&mut self) -> Result<Option<(Line, Line, Line)>, FastqError> {
        self.rec.clear();
        // Repacking (buffer shift) is allowed only here, at the record
        // boundary: no earlier line spans are in use. Every later line
        // read must not shift the buffer, because the header/sequence
        // offsets are still needed for the deferred checks below.
        let Some(hdr) = self.read_line(true)? else {
            return Ok(None);
        };
        self.ascii_check(&hdr)?;
        let Some(seq) = self.read_line(false)? else {
            return Err(self.premature_eof());
        };
        self.ascii_check(&seq)?;
        let hfirst_ok =
            hdr.start < hdr.end && (self.buf[hdr.start] == b'@' || self.buf[hdr.start] == b'>');
        if !hfirst_ok {
            return Err(format_error(
                Some(hdr.no),
                format!(
                    "Line expected to start with '@', but found '{}'",
                    Self::found_char(&hdr, &self.buf)
                ),
            ));
        }
        let Some(plus) = self.read_line(false)? else {
            return Err(self.premature_eof());
        };
        self.ascii_check(&plus)?;
        if plus.start == plus.end || self.buf[plus.start] != b'+' {
            return Err(format_error(
                Some(plus.no),
                format!(
                    "Line expected to start with '+', but found '{}'",
                    Self::found_char(&plus, &self.buf)
                ),
            ));
        }
        let Some(qual) = self.read_line(false)? else {
            return Err(self.premature_eof());
        };
        self.ascii_check(&qual)?;
        if (seq.end - seq.start) != (qual.end - qual.start) {
            return Err(format_error(
                Some(qual.no),
                "Length of sequence and qualities differ".to_string(),
            ));
        }
        Ok(Some((hdr, seq, qual)))
    }
}

impl<R: Read> FastqReader<R> {
    /// Consume one record, returning the name0 slice (header minus the
    /// leading `@`/`>`, up to the first space byte). The slice is valid
    /// until the next call.
    pub fn next_name0(&mut self) -> Option<Result<&[u8], FastqError>> {
        Some(match self.parse_record() {
            Ok(None) => return None,
            Ok(Some((hdr, _seq, _qual))) => {
                let name = &self.buf[(hdr.start + 1)..hdr.end];
                let n0 = name.iter().position(|b| *b == b' ').unwrap_or(name.len());
                Ok(&self.buf[(hdr.start + 1)..(hdr.start + 1 + n0)])
            }
            Err(e) => Err(e),
        })
    }

    /// Consume one record, returning the owned normalised record.
    pub fn next_record(&mut self) -> Option<Result<Record, FastqError>> {
        Some(match self.parse_record() {
            Ok(None) => return None,
            Ok(Some((hdr, seq, qual))) => Ok(Record {
                name: self.buf[(hdr.start + 1)..hdr.end].to_vec(),
                seq: self.buf[seq.start..seq.end].to_vec(),
                qual: self.buf[qual.start..qual.end].to_vec(),
            }),
            Err(e) => Err(e),
        })
    }
}

/// A file-backed reader: gzip (multi-member) when the extension is `.gz`,
/// plain otherwise — the same extension-based dispatch as Python's xopen.
pub enum FileFastqReader {
    Plain(FastqReader<File>),
    /// `MultiGzDecoder`: reads concatenated gzip members, like Python's
    /// gzip module (a plain `GzDecoder` stops at the first member).
    Gz(FastqReader<MultiGzDecoder<File>>),
}

impl FileFastqReader {
    pub fn open(path: &Path) -> io::Result<Self> {
        let is_gz = path
            .extension()
            .and_then(|e| e.to_str())
            .is_some_and(|e| e == "gz");
        if is_gz {
            Ok(Self::Gz(FastqReader::new(MultiGzDecoder::new(File::open(
                path,
            )?))))
        } else {
            Ok(Self::Plain(FastqReader::new(File::open(path)?)))
        }
    }

    pub fn next_name0(&mut self) -> Option<Result<&[u8], FastqError>> {
        match self {
            Self::Plain(r) => r.next_name0(),
            Self::Gz(r) => r.next_name0(),
        }
    }

    pub fn next_record(&mut self) -> Option<Result<Record, FastqError>> {
        match self {
            Self::Plain(r) => r.next_record(),
            Self::Gz(r) => r.next_record(),
        }
    }
}
