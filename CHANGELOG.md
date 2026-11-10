## Unreleased

### Feat

- Rust rewrite of the engine (`fastq-purge-rs`, `rust/fastq-purge/`; all
  milestones merged into master): canonical semantics per `RUST_REWRITE.md` —
  both methods purge the undetermined files, comparing by name0 (D1, D2);
  new `--purge-assigned` flag inverts the direction (D3); the approx method
  is rebuilt on a memory-bounded per-lane bloom filter with optimal probe
  count (D4), whose purged outputs are byte-identical to the exact method
  for any `--fpr`. Pipelines are threaded (rayon) and differentially
  verified against the Python implementation (byte-identical purged
  outputs). Side-by-side rollout, verified 2026-11-10 on a real NovaSeqX
  flowcell: purged outputs byte-identical after decompression, logs differ
  only in documented ways (`RUST_REWRITE.md` §8). The Python CLI will be
  deprecated once the rollout completes.

### Fix

- `--method approx` CLI path: repaired the `AttributeError` crash (a
  lane→pair dict was passed where the bloom builder expects the undetermined
  R1 file list, and the output directory was lost to the per-lane output
  dict). The bloom hash is now the builtin `hash` instead of
  `sha256(pickle(...))` — filter build is ~3x faster, since build and query
  happen in one process. Semantics unchanged: approx still purges the
  assigned files keyed by the full record name; the exact/approx divergence
  is documented, not changed, in the Python package (`RUST_REWRITE.md` §2).

### Note

- The legacy mode (`--undetermined-path`/`--assigned-path`) accumulates the
  undetermined set across lanes (module global at `main.py:37`;
  `build_undetermined_set` mutates it and each per-lane call rebinds to it),
  so the per-lane `Number of undetermined reads` and `Found {N} duplicates in
  the undetermined file.` lines are cumulative across lanes (e.g.
  10000/20000/30000/40000 over a 4-lane run). Purged outputs are unaffected:
  a record removed from a lane's file carries a name from that file, hence
  already in that lane's own set, so the cumulative set changes the reported
  counts, not the purged records. Not fixed in this package.
  `fastq-purge-rs` uses fresh per-lane sets (correct per-lane counts,
  identical outputs). The `samplesheet_support` branch has the same pattern
  in sheet mode (D9); the legacy-mode instance is D10
  (`RUST_REWRITE.md` §3).
- The legacy mode has no `Undetermined` file exclusion: pointing
  `--assigned-path` at a directory that also contains the undetermined files
  makes every undetermined read match itself (100% purge, empty outputs).
  Pass only the demultiplexed (assigned) files; the sheet mode excludes
  `Undetermined` files from the assigned set.

## 0.2.0 (2025-05-23)

### Feat

- Add `commitizen` section to `pyproject.toml`
- Add `pre-commit` and `commitizen` as `dev` dependencies
- Introduce semaphore to limit process submissions
- Implement parallelisation
- Add option for tesing `loky` and `multiprocess` side by side
- Switch to `dnaio` for reading and writing
- Include `psutil` for memory usage tracking
- First working version of the code
- Initial commit of main repository files

### Fix

- Arguments and program name fix

### Refactor

- Large refactoring, collecting code into functions and classes
- Extensive refactoring of the code
