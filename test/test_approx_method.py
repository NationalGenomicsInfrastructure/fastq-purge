"""
Regression tests for the ``--method approx`` CLI path.

The released CLI used to crash with ``AttributeError: 'str' object has no
attribute 'name'`` because ``main()`` passed the per-lane argument dictionaries
straight into ``build_bloom_filter`` / ``process_target_file``, which expect
lists of file paths.
"""

import re
import subprocess
import sys

import dnaio


def write_fastq(path, names, length=150):
    """Write a gzipped FASTQ file whose records carry the given full names."""
    with dnaio.open(path, mode="w", compression_level=1) as writer:
        for name in names:
            record = dnaio.SequenceRecord(name, "A" * length, None)
            record.qualities = "I" * length
            writer.write(record)


def undetermined_names(n):
    return [
        f"M00150:2:AABCCC:1:1101:{i:04d}:{i:04d} 1:N:0:UND{i:04d}+UND{i:04d}"
        for i in range(n)
    ]


def sample_names(n, dup_names, n_dup):
    """First ``n_dup`` names are byte-identical to ``dup_names``; the rest unique."""
    names = list(dup_names[:n_dup])
    for i in range(n_dup, n):
        names.append(
            f"M00150:2:AABCCC:1:1101:{i:04d}:{i:04d} 1:N:0:SMP{i:04d}+SMP{i:04d}"
        )
    return names


ANSI_RE = re.compile(r"\x1b\[[0-9;]*m")


def test_approx_method_cli(tmp_path):
    n_und, n_total, n_dup = 100, 50, 30

    und_dir = tmp_path / "und"
    und_dir.mkdir()
    smp_dir = tmp_path / "samples"
    smp_dir.mkdir()
    out_dir = tmp_path / "out"

    und_names = undetermined_names(n_und)
    # Paired reads share their full names
    write_fastq(und_dir / "Undetermined_S0_L001_R1_001.fastq.gz", und_names)
    write_fastq(und_dir / "Undetermined_S0_L001_R2_001.fastq.gz", und_names)

    smp_names = sample_names(n_total, und_names, n_dup)
    write_fastq(smp_dir / "S1_L001_R1_001.fastq.gz", smp_names)
    write_fastq(smp_dir / "S1_L001_R2_001.fastq.gz", smp_names)

    proc = subprocess.run(
        [
            sys.executable,
            "-m",
            "fastq_purge.main",
            "--undetermined-path",
            str(und_dir),
            "--assigned-path",
            str(smp_dir),
            "--output-path",
            str(out_dir),
            "--method",
            "approx",
        ],
        capture_output=True,
        text=True,
    )
    assert proc.returncode == 0, proc.stderr

    stdout = ANSI_RE.sub("", proc.stdout)
    assert f"Removed {n_dup} reads" in stdout
    assert f"Kept {n_total - n_dup} reads" in stdout

    purged = out_dir / "S1_L001_R1_001.purged.fastq.gz"
    assert purged.is_file()
    kept = [rec.name for rec in dnaio.open(purged)]
    assert len(kept) == n_total - n_dup
    # Only the reads absent from the undetermined data may survive
    assert set(kept).isdisjoint(set(und_names))
    assert set(kept) == set(smp_names[n_dup:])
