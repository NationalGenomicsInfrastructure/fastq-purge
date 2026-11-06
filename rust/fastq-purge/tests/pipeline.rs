//! End-to-end tests for the exact pipeline (milestone M3): the binary is
//! driven with tiny synthetic datasets and the purged outputs + log
//! lines are checked against the canonical contract (RUST_REWRITE.md §3,
//! decisions D1–D3).

use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> &'static str {
    env!("CARGO_BIN_EXE_fastq-purge-rs")
}

fn tmpdir(tag: &str) -> PathBuf {
    let d = std::env::temp_dir().join(format!("fqp_m3_{}_{}", tag, std::process::id()));
    let _ = std::fs::remove_dir_all(&d);
    std::fs::create_dir_all(&d).unwrap();
    d
}

/// One record per name: `@{name} 1:N:0\nACGT\n+\nIIII\n`.
fn write_fastq(path: &Path, names: &[&str]) {
    let mut s = String::new();
    for n in names {
        s.push_str(&format!("@{n} 1:N:0\nACGT\n+\nIIII\n"));
    }
    std::fs::write(path, s).unwrap();
}

/// Run the binary with the given extra flags from directory `d`.
fn run(d: &Path, extra: &[&str]) -> String {
    let mut args: Vec<&str> = vec![
        "--undetermined-path",
        "und",
        "--assigned-path",
        "assigned",
        "--output-path",
        "out",
    ];
    args.extend_from_slice(extra);
    let out = Command::new(bin())
        .current_dir(d)
        .args(&args)
        .output()
        .unwrap();
    let stdout = String::from_utf8_lossy(&out.stdout).to_string();
    assert_eq!(
        out.status.code().unwrap_or(-1),
        0,
        "stdout:\n{stdout}\nstderr:\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    stdout
}

/// Standard 4-read undetermined lane (L001) + one assigned file.
fn setup(d: &Path, assigned: &[&str], r2: bool) {
    std::fs::create_dir_all(d.join("und")).unwrap();
    std::fs::create_dir_all(d.join("assigned")).unwrap();
    std::fs::create_dir_all(d.join("out")).unwrap();
    write_fastq(
        &d.join("und/Undetermined_S0_L001_R1_001.fastq"),
        &["a", "b", "c", "d"],
    );
    if r2 {
        write_fastq(
            &d.join("und/Undetermined_S0_L001_R2_001.fastq"),
            &["e", "f", "g", "h"],
        );
    }
    write_fastq(&d.join("assigned/S1_L001_R1_001.fastq"), assigned);
}

#[test]
fn paired_lane_purges_both_outputs() {
    let d = tmpdir("paired");
    setup(&d, &["b", "x", "d"], true);
    let stdout = run(&d, &[]);
    assert!(stdout.contains("Found 2 duplicates in the undetermined file."));
    assert!(stdout.contains(
        "Writing purged fastq files 'Undetermined_S0_L001_R1_001.purged.fastq' and 'Undetermined_S0_L001_R2_001.purged.fastq'"
    ));
    assert!(stdout.contains("    Removed 2 reads"));
    assert!(stdout.contains("    Kept 2 reads"));
    assert!(stdout.contains("Done!"));

    let p1 = d.join("out/Undetermined_S0_L001_R1_001.purged.fastq");
    let p2 = d.join("out/Undetermined_S0_L001_R2_001.purged.fastq");
    assert_eq!(
        std::fs::read_to_string(&p1).unwrap(),
        "@a 1:N:0\nACGT\n+\nIIII\n@c 1:N:0\nACGT\n+\nIIII\n"
    );
    // paired: R2 keeps the records whose R1 name0 was kept
    assert_eq!(
        std::fs::read_to_string(&p2).unwrap(),
        "@e 1:N:0\nACGT\n+\nIIII\n@g 1:N:0\nACGT\n+\nIIII\n"
    );
}

#[test]
fn no_duplicates_writes_nothing() {
    let d = tmpdir("nodups");
    setup(&d, &["z", "y"], true);
    let stdout = run(&d, &[]);
    assert!(stdout.contains("Found 0 duplicates in the undetermined file."));
    assert!(stdout.contains("Found no duplicates in the undetermined reads. Exiting..."));
    assert!(stdout.contains("Done!"));
    assert!(
        std::fs::read_dir(d.join("out")).unwrap().count() == 0,
        "no output files may be created when nothing is purged"
    );
}

#[test]
fn purge_assigned_flag_inverts_direction() {
    let d = tmpdir("purgeassigned");
    setup(&d, &["b", "x", "d"], true);
    let stdout = run(&d, &["--purge-assigned"]);
    assert!(stdout.contains("Purging assigned files..."));
    assert!(stdout.contains("    Removed 2 reads"));
    assert!(stdout.contains("    Kept 1 reads"));
    let pa = d.join("out/S1_L001_R1_001.purged.fastq");
    assert_eq!(
        std::fs::read_to_string(&pa).unwrap(),
        "@x 1:N:0\nACGT\n+\nIIII\n"
    );
    // no undetermined purged outputs in this mode
    assert!(!d
        .join("out/Undetermined_S0_L001_R1_001.purged.fastq")
        .exists());
}

#[test]
fn single_read_lane() {
    let d = tmpdir("single");
    setup(&d, &["b", "d"], false);
    let stdout = run(&d, &[]);
    assert!(stdout.contains("Writing purged fastq file 'Undetermined_S0_L001_R1_001.purged.fastq'"));
    assert!(stdout.contains("    Removed 2 reads"));
    assert!(stdout.contains("    Kept 2 reads"));
    let p1 = d.join("out/Undetermined_S0_L001_R1_001.purged.fastq");
    assert_eq!(
        std::fs::read_to_string(&p1).unwrap(),
        "@a 1:N:0\nACGT\n+\nIIII\n@c 1:N:0\nACGT\n+\nIIII\n"
    );
    assert!(!d
        .join("out/Undetermined_S0_L001_R2_001.purged.fastq")
        .exists());
}

#[test]
fn gzipped_end_to_end() {
    // the production path: gzipped inputs, gzipped purged outputs
    let d = tmpdir("gz");
    std::fs::create_dir_all(d.join("und")).unwrap();
    std::fs::create_dir_all(d.join("assigned")).unwrap();
    std::fs::create_dir_all(d.join("out")).unwrap();
    for (dir, name, names) in [
        (
            "und",
            "Undetermined_S0_L001_R1_001.fastq.gz",
            &["a", "b", "c"][..],
        ),
        ("assigned", "S1_L001_R1_001.fastq.gz", &["b", "x"][..]),
    ] {
        let mut s = String::new();
        for n in names {
            s.push_str(&format!("@{n} 1:N:0\nACGT\n+\nIIII\n"));
        }
        let mut enc = flate2::write::GzEncoder::new(
            std::fs::File::create(d.join(dir).join(name)).unwrap(),
            flate2::Compression::default(),
        );
        std::io::Write::write_all(&mut enc, s.as_bytes()).unwrap();
        enc.finish().unwrap();
    }
    let stdout = run(&d, &[]);
    assert!(stdout.contains("Found 1 duplicates in the undetermined file."));
    let out = d.join("out/Undetermined_S0_L001_R1_001.purged.fastq.gz");
    let mut dec = flate2::read::GzDecoder::new(std::fs::File::open(&out).unwrap());
    let mut got = String::new();
    std::io::Read::read_to_string(&mut dec, &mut got).unwrap();
    assert_eq!(got, "@a 1:N:0\nACGT\n+\nIIII\n@c 1:N:0\nACGT\n+\nIIII\n");
}

// ---- milestone M4: the approx (bloom) pipeline -------------------------

/// The D4 guarantee (RUST_REWRITE.md §8/M4): for any `--fpr` the purged
/// outputs are identical to the exact method's; only the `Found {K}`
/// count may differ (by the number of structural false positives).
#[test]
fn approx_matches_exact_output() {
    let d = tmpdir("approx_exact");
    setup(&d, &["b", "x", "d"], true);
    let exact_out = run(&d, &[]);

    // snapshot the exact purged outputs, then clear the out dir
    let p1 = d.join("out/Undetermined_S0_L001_R1_001.purged.fastq");
    let p2 = d.join("out/Undetermined_S0_L001_R2_001.purged.fastq");
    let exact_files = (std::fs::read(&p1).unwrap(), std::fs::read(&p2).unwrap());
    std::fs::remove_dir_all(d.join("out")).unwrap();
    std::fs::create_dir_all(d.join("out")).unwrap();

    let approx_out = run(
        &d,
        &[
            "--method",
            "approx",
            "--max-items",
            "100000",
            "--fpr",
            "0.01",
        ],
    );
    assert!(approx_out.contains("Using bloom filter method"));
    assert!(approx_out.contains("Building bloom filter..."));
    assert!(approx_out.contains("Scanning assigned files for already-assigned reads..."));

    // same Removed/Kept counts as the exact run
    for line in ["    Removed 2 reads", "    Kept 2 reads", "Done!"] {
        assert!(approx_out.contains(line), "missing {line:?}\n{approx_out}");
        assert!(exact_out.contains(line), "missing {line:?}\n{exact_out}");
    }

    // byte-identical purged outputs (the D4 guarantee)
    assert_eq!(std::fs::read(&p1).unwrap(), exact_files.0);
    assert_eq!(std::fs::read(&p2).unwrap(), exact_files.1);
}

#[test]
fn approx_no_duplicates_writes_nothing() {
    let d = tmpdir("approx_nodups");
    setup(&d, &["z", "y"], true);
    let stdout = run(
        &d,
        &[
            "--method",
            "approx",
            "--max-items",
            "100000",
            "--fpr",
            "0.01",
        ],
    );
    assert!(stdout.contains("Found 0 duplicates in the undetermined file."));
    assert!(stdout.contains("Found no duplicates in the undetermined reads. Exiting..."));
    assert!(stdout.contains("Done!"));
    assert_eq!(std::fs::read_dir(d.join("out")).unwrap().count(), 0);
}

#[test]
fn approx_purge_assigned_flag() {
    let d = tmpdir("approx_purgeassigned");
    setup(&d, &["b", "x", "d"], true);
    let stdout = run(
        &d,
        &[
            "--method",
            "approx",
            "--max-items",
            "100000",
            "--fpr",
            "0.01",
            "--purge-assigned",
        ],
    );
    assert!(stdout.contains("Purging assigned files..."));
    assert!(stdout.contains("    Removed 2 reads"));
    assert!(stdout.contains("    Kept 1 reads"));
    let pa = d.join("out/S1_L001_R1_001.purged.fastq");
    assert_eq!(
        std::fs::read_to_string(&pa).unwrap(),
        "@x 1:N:0\nACGT\n+\nIIII\n"
    );
    assert!(!d
        .join("out/Undetermined_S0_L001_R1_001.purged.fastq")
        .exists());
}
