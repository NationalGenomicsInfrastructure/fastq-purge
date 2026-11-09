//! Integration tests for the FASTQ reader/writer: the full empirical
//! matrix from the dnaio 1.2.3 study (RUST_REWRITE.md §8, M2 findings).
//! Every expected error string and line number was observed from the
//! Python implementation.

use std::io::{Cursor, Write};

use fastq_purge::fastq::{
    FastqError, FastqReader, FastqWriter, FileFastqReader, PairedFastqWriter, Record,
};

fn records(data: &[u8]) -> Vec<Result<Record, FastqError>> {
    let mut r = FastqReader::new(Cursor::new(data.to_vec()));
    let mut out = Vec::new();
    while let Some(x) = r.next_record() {
        out.push(x);
    }
    out
}

fn names(data: &[u8]) -> Result<Vec<Vec<u8>>, FastqError> {
    let mut r = FastqReader::new(Cursor::new(data.to_vec()));
    let mut out = Vec::new();
    while let Some(x) = r.next_name0() {
        out.push(x?.to_vec());
    }
    Ok(out)
}

fn first_err(data: &[u8]) -> String {
    let mut r = FastqReader::new(Cursor::new(data.to_vec()));
    loop {
        match r.next_record() {
            Some(Ok(_)) => {}
            Some(Err(e)) => return e.to_string(),
            None => panic!("expected an error, got clean EOF"),
        }
    }
}

fn tmpdir() -> std::path::PathBuf {
    let d = std::env::temp_dir().join(format!("fqp_m2_{}", std::process::id()));
    std::fs::create_dir_all(&d).unwrap();
    d
}

// ---------------------------------------------------------------------------
// happy path + normalisation
// ---------------------------------------------------------------------------

#[test]
fn valid_basic_record() {
    let recs = records(b"@read1 extra stuff\nACGT\n+\nIIII\n");
    assert_eq!(recs.len(), 1);
    let r = recs[0].clone().unwrap();
    assert_eq!(r.name, b"read1 extra stuff");
    assert_eq!(r.seq, b"ACGT");
    assert_eq!(r.qual, b"IIII");
    assert_eq!(r.name0(), b"read1");
    assert_eq!(r.normalized_bytes(), b"@read1 extra stuff\nACGT\n+\nIIII\n");
}

#[test]
fn crlf_is_normalised() {
    let recs = records(b"@read2 x\r\nACGT\r\n+\r\nIIII\r\n");
    let r = recs[0].clone().unwrap();
    assert_eq!(r.name, b"read2 x");
    assert_eq!(r.seq, b"ACGT");
    assert_eq!(r.qual, b"IIII");
    assert_eq!(r.normalized_bytes(), b"@read2 x\nACGT\n+\nIIII\n");
}

#[test]
fn named_plus_line_is_dropped() {
    let recs = records(b"@read1 extra stuff\nACGT\n+read1 extra stuff\nIIII\n");
    let r = recs[0].clone().unwrap();
    assert_eq!(r.normalized_bytes(), b"@read1 extra stuff\nACGT\n+\nIIII\n");
}

#[test]
fn name0_stops_at_first_space() {
    assert_eq!(
        names(b"@a 1:N:0\nACGT\n+\nIIII\n").unwrap(),
        vec![b"a".to_vec()]
    );
}

#[test]
fn tab_is_not_a_separator() {
    assert_eq!(
        names(b"@read12\tACGT\nACGT\n+\nIIII\n").unwrap(),
        vec![b"read12\tACGT".to_vec()]
    );
}

#[test]
fn no_space_name0_is_whole_name() {
    assert_eq!(
        names(b"@nospacetoken\nACGT\n+\nIIII\n").unwrap(),
        vec![b"nospacetoken".to_vec()]
    );
}

#[test]
fn empty_name() {
    let recs = records(b"@\nACGT\n+\nIIII\n");
    let r = recs[0].clone().unwrap();
    assert_eq!(r.name, b"");
    assert_eq!(r.name0(), b"");
}

#[test]
fn no_trailing_newline_last_record() {
    let recs = records(b"@a x\nACGT\n+\nIIII");
    assert_eq!(recs.len(), 1);
    assert!(recs[0].is_ok());
}

#[test]
fn empty_file_is_zero_records() {
    assert!(records(b"").is_empty());
}

#[test]
fn two_records() {
    let data = b"@a x\nACGT\n+\nIIII\n@b y\nTTTT\n+\nJJJJ\n";
    let recs = records(data);
    assert_eq!(recs.len(), 2);
    assert_eq!(names(data).unwrap(), vec![b"a".to_vec(), b"b".to_vec()]);
}

#[test]
fn gt_header_is_accepted_like_at() {
    // Documented deviation: dnaio 1.2.3 mis-parses `>`-header files
    // (two-header misdetection); here they are treated as ordinary records.
    let recs = records(b">read3\nACGT\n+\nIIII\n");
    assert_eq!(recs.len(), 1);
    let r = recs[0].clone().unwrap();
    assert_eq!(r.name, b"read3");
    assert_eq!(r.seq, b"ACGT");
}

// ---------------------------------------------------------------------------
// error strings and line numbers (exact dnaio parity)
// ---------------------------------------------------------------------------

#[test]
fn truncated_three_lines() {
    assert_eq!(
        first_err(b"@a x\nACGT\n+\n"),
        "Error in FASTQ file at line 4: Premature end of file \
         encountered. The incomplete final record was: '@a x\\nACGT\\n+\\n'"
    );
}

#[test]
fn partial_record_at_end() {
    assert_eq!(
        first_err(b"@a x\nACGT\n+\nIIII\n@b y\n"),
        "Error in FASTQ file at line 6: Premature end of file \
         encountered. The incomplete final record was: '@b y\\n'"
    );
}

#[test]
fn trailing_blank_line() {
    assert_eq!(
        first_err(b"@a x\nACGT\n+\nIIII\n\n"),
        "Error in FASTQ file at line 6: Premature end of file \
         encountered. The incomplete final record was: '\\n'"
    );
}

#[test]
fn missing_at_second_record() {
    assert_eq!(
        first_err(b"@a x\nACGT\n+\nIIII\nb noat\nACGT\n+\nIIII\n"),
        "Error in FASTQ file at line 5: Line expected to start with '@', \
         but found 'b'"
    );
}

#[test]
fn blank_line_between_records() {
    assert_eq!(
        first_err(b"@a x\nACGT\n+\nIIII\n\n@b y\nTTTT\n+\nIIII\n"),
        "Error in FASTQ file at line 5: Line expected to start with '@', \
         but found '\\n'"
    );
}

#[test]
fn ws_only_line_between_records() {
    assert_eq!(
        first_err(b"@a x\nACGT\n+\nIIII\n   \n@b y\nTTTT\n+\nIIII\n"),
        "Error in FASTQ file at line 5: Line expected to start with '@', \
         but found ' '"
    );
}

#[test]
fn plus_line_must_start_with_plus() {
    assert_eq!(
        first_err(b"@a x\nACGT\nxIIII\nIIII\n"),
        "Error in FASTQ file at line 3: Line expected to start with '+', \
         but found 'x'"
    );
}

#[test]
fn seq_qual_length_checks() {
    let msg = "Error in FASTQ file at line 4: Length of sequence and \
               qualities differ";
    assert_eq!(first_err(b"@a x\n\n+\nIIII\n"), msg); // empty seq
    assert_eq!(first_err(b"@a x\nACGT\n+\n\n"), msg); // empty qual
    assert_eq!(first_err(b"@a x\nACGT\n+\nII\n"), msg); // qual short
    assert_eq!(first_err(b"@a x\nAC\n+\nIIIIII\n"), msg); // qual long
    assert_eq!(first_err(b"@a x\nAC GT\n+\nIIII\n"), msg); // space in seq
}

#[test]
fn non_ascii_any_line() {
    let msg = "Error in FASTQ file at unknown line: Non-ASCII characters \
               found in record.";
    assert_eq!(first_err(b"@\xe9bad x\nACGT\n+\nIIII\n"), msg); // name
    assert_eq!(first_err(b"@a x\n\xe9CGT\n+\nIIII\n"), msg); // seq
    assert_eq!(first_err(b"@a x\nACGT\n+\xe9\nIIII\n"), msg); // plus content
    assert_eq!(first_err(b"@a x\nACGT\n+\nII\xe9I\n"), msg); // qual
}

#[test]
fn non_ascii_fires_before_truncation() {
    let msg = "Error in FASTQ file at unknown line: Non-ASCII characters \
               found in record.";
    assert_eq!(first_err(b"@\xe9bad x\nACGT\n+\n"), msg);
}

#[test]
fn incomplete_record_quote_escaping() {
    // single quote in content -> double-quoted Python repr
    assert_eq!(
        first_err(b"@it's x\nACGT\n+\n"),
        "Error in FASTQ file at line 4: Premature end of file encountered. \
         The incomplete final record was: \"@it's x\\nACGT\\n+\\n\""
    );
    // double quote in content -> stays single-quoted
    assert_eq!(
        first_err(b"@a\"b x\nACGT\n+\n"),
        "Error in FASTQ file at line 4: Premature end of file encountered. \
         The incomplete final record was: '@a\"b x\\nACGT\\n+\\n'"
    );
    // CRLF raw bytes quoted with escapes; tab escaped
    assert_eq!(
        first_err(b"@a x\r\nACGT\r\n+\r\n"),
        "Error in FASTQ file at line 4: Premature end of file encountered. \
         The incomplete final record was: '@a x\\r\\nACGT\\r\\n+\\r\\n'"
    );
    assert_eq!(
        first_err(b"@a\tb\nACGT\n+\n"),
        "Error in FASTQ file at line 4: Premature end of file encountered. \
         The incomplete final record was: '@a\\tb\\nACGT\\n+\\n'"
    );
}

// ---------------------------------------------------------------------------
// file-backed reading: plain, gzip, multi-member gzip
// ---------------------------------------------------------------------------

#[test]
fn open_plain_file() {
    let d = tmpdir();
    let p = d.join("p.fastq");
    std::fs::write(&p, b"@a x\nACGT\n+\nIIII\n@b y\nTTTT\n+\nJJJJ\n").unwrap();
    let mut r = FileFastqReader::open(&p).unwrap();
    assert_eq!(r.next_name0().unwrap().unwrap(), b"a");
    assert_eq!(r.next_name0().unwrap().unwrap(), b"b");
    assert!(r.next_name0().is_none());
}

#[test]
fn open_gz_file() {
    let d = tmpdir();
    let p = d.join("g.fastq.gz");
    let mut enc = flate2::write::GzEncoder::new(
        std::fs::File::create(&p).unwrap(),
        flate2::Compression::default(),
    );
    enc.write_all(b"@a x\nACGT\n+\nIIII\n").unwrap();
    enc.finish().unwrap();
    let mut r = FileFastqReader::open(&p).unwrap();
    assert_eq!(r.next_name0().unwrap().unwrap(), b"a");
    assert!(r.next_name0().is_none());
}

#[test]
fn multi_member_gz_is_supported() {
    // Python's gzip module reads concatenated members; the reader must too.
    let d = tmpdir();
    let p = d.join("multi.fastq.gz");
    let f = std::fs::File::create(&p).unwrap();
    let bw = std::io::BufWriter::new(f);
    let mut enc1 = flate2::write::GzEncoder::new(bw, flate2::Compression::default());
    enc1.write_all(b"@a x\nACGT\n+\nIIII\n").unwrap();
    let bw = enc1.finish().unwrap();
    let mut enc2 = flate2::write::GzEncoder::new(bw, flate2::Compression::default());
    enc2.write_all(b"@b y\nTTTT\n+\nIIII\n").unwrap();
    enc2.finish().unwrap().into_inner().unwrap();

    let mut r = FileFastqReader::open(&p).unwrap();
    assert_eq!(r.next_name0().unwrap().unwrap(), b"a");
    assert_eq!(r.next_name0().unwrap().unwrap(), b"b");
    assert!(r.next_name0().is_none());
}

// ---------------------------------------------------------------------------
// writers
// ---------------------------------------------------------------------------

#[test]
fn plain_writer_exact_bytes() {
    let d = tmpdir();
    let p = d.join("w.fastq");
    let mut w = FastqWriter::create(&p).unwrap();
    w.write_record(b"z w", b"ACGT", b"IIII").unwrap();
    w.write_record_bytes(b"@q r\nTTTT\n+\nJJJJ\n").unwrap();
    w.finish().unwrap();
    assert_eq!(
        std::fs::read(&p).unwrap(),
        b"@z w\nACGT\n+\nIIII\n@q r\nTTTT\n+\nJJJJ\n"
    );
}

#[test]
fn gz_writer_decompressed_bytes() {
    let d = tmpdir();
    let p = d.join("w.fastq.gz");
    let mut w = FastqWriter::create(&p).unwrap();
    w.write_record(b"z w", b"ACGT", b"IIII").unwrap();
    w.finish().unwrap();
    let raw = std::fs::read(&p).unwrap();
    // gzip magic + (xopen-style) zero mtime
    assert_eq!(&raw[0..2], b"\x1f\x8b");
    let mut dec = flate2::read::GzDecoder::new(&raw[..]);
    let mut out = Vec::new();
    std::io::Read::read_to_end(&mut dec, &mut out).unwrap();
    assert_eq!(out, b"@z w\nACGT\n+\nIIII\n");
}

#[test]
fn paired_writer_lockstep() {
    let d = tmpdir();
    let p1 = d.join("r1.fastq");
    let p2 = d.join("r2.fastq");
    let mut w = PairedFastqWriter::create(&p1, &p2).unwrap();
    w.write_pair(
        &Record {
            name: b"p1 x".to_vec(),
            seq: b"ACGT".to_vec(),
            qual: b"IIII".to_vec(),
        },
        &Record {
            name: b"p2 x".to_vec(),
            seq: b"TTTT".to_vec(),
            qual: b"JJJJ".to_vec(),
        },
    )
    .unwrap();
    w.finish().unwrap();
    assert_eq!(std::fs::read(&p1).unwrap(), b"@p1 x\nACGT\n+\nIIII\n");
    assert_eq!(std::fs::read(&p2).unwrap(), b"@p2 x\nTTTT\n+\nJJJJ\n");
}

#[test]
fn writer_round_trip_normalisation() {
    // read a CRLF, named-plus file; write it back; bytes are normalised
    let d = tmpdir();
    let src = d.join("in.fastq");
    std::fs::write(
        &src,
        b"@read1 extra stuff\r\nACGT\r\n+read1 extra stuff\r\nIIII\r\n",
    )
    .unwrap();
    let dst = d.join("out.fastq");
    {
        let mut r = FileFastqReader::open(&src).unwrap();
        let mut w = FastqWriter::create(&dst).unwrap();
        while let Some(rec) = r.next_record() {
            let rec = rec.unwrap();
            w.write_record_bytes(&rec.normalized_bytes()).unwrap();
        }
        w.finish().unwrap();
    }
    assert_eq!(
        std::fs::read(&dst).unwrap(),
        b"@read1 extra stuff\nACGT\n+\nIIII\n"
    );
}
