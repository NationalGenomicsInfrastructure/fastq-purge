## Unreleased

### Feat

- Rust rewrite of the engine (`fastq-purge-rs`, `rust/fastq-purge/`; branch
  `rust/m1-cli-discovery`): canonical semantics per `RUST_REWRITE.md` — both
  methods purge the undetermined files, comparing by name0 (D1, D2); new
  `--purge-assigned` flag inverts the direction (D3); the approx method is
  rebuilt on a memory-bounded per-lane bloom filter with optimal probe count
  (D4), whose purged outputs are byte-identical to the exact method for any
  `--fpr`. Pipelines are threaded (rayon) and differentially verified against
  the Python implementation (byte-identical purged outputs). Ships
  side-by-side; see `RUST_REWRITE.md` §8 for the rollout plan.

### Fix

- `--method approx` CLI path: repaired the `AttributeError` crash (a
  lane→pair dict was passed where the bloom builder expects the undetermined
  R1 file list, and the output directory was lost to the per-lane output
  dict). The bloom hash is now the builtin `hash` instead of
  `sha256(pickle(...))` — filter build is ~3x faster, since build and query
  happen in one process. Semantics unchanged: approx still purges the
  assigned files keyed by the full record name; the exact/approx divergence
  is documented, not changed, in the Python package (`RUST_REWRITE.md` §2).

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
