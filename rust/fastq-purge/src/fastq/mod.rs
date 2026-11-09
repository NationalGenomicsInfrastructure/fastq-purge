//! Streaming FASTQ(.gz) reader and writer.
//!
//! The reader is a line-based port of dnaio 1.2.3's C parser, verified
//! empirically against it (see RUST_REWRITE.md §8, M2 findings):
//!
//! - a record is strictly 4 lines: header (`@` or `>`), sequence, `+`, qual
//! - CRLF is normalised to LF; the `+` line's content (a repeated name) is
//!   dropped — the writer always emits a bare `+`
//! - `name0` is the header (minus the leading `@`/`>`) up to the first
//!   space byte (0x20); tabs are not separators; no space => whole name
//! - error model (exact messages, 1-based physical line numbers):
//!   * header check is DEFERRED until the sequence line has been read: a
//!     bad header followed by EOF reports premature-EOF, not the header
//!     error (matches the C parser's behaviour)
//!   * the `+` check fires when the `+` line is read
//!   * `Length of sequence and qualities differ` fires after the qual line,
//!     reported at the qual line number
//!   * non-ASCII bytes abort as soon as the offending line is read
//!     ("Error in FASTQ file at unknown line: Non-ASCII characters found
//!     in record.")
//!   * premature-EOF quotes the raw bytes of the incomplete record attempt
//!     in Python-repr style
//! - gzip input: multi-member aware (like Python's gzip module)

pub mod read;
pub mod write;

pub use read::{FastqError, FastqReader, FileFastqReader, Record};
pub use write::{FastqWriter, PairedFastqWriter};
