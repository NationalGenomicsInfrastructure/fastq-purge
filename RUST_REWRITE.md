# fastq-purge — Rust rewrite design

Status: DRAFT for review · Profiled: 2026-11-04 · Scope: full rewrite of the `fastq-purge` package as a standalone Rust binary.

## 1. Motivation & measured baseline

Current state (Python 3.11, `src/fastq_purge/main.py`, 869 lines, single file):

| Phase (5M-read synthetic set, 20% engineered duplicates) | Time | Per record |
|---|---|---|
| exact: build set from undetermined R1 (`build_undetermined_set`) | 2.93 s | 0.59 µs |
| exact: scan assigned files (`pool_task`) | 3.96 s | 0.79 µs |
| exact: re-read undetermined, write purged (`write_purged_fastq`) | 2.21 s | 0.44 µs |
| **exact total (single process)** | **9.1 s** | ~1.8 µs |
| approx: build bloom from undetermined (`build_bloom_filter`) | 7.96 s | 1.59 µs |
| approx: filter + write (`process_target_file`, gzip L7) | 30.46 s | 6.09 µs |
| **approx total** | **38.4 s** | ~7.7 µs |

Memory: Python set of name-strings ≈ 90–120 B/name → **~9–12 GB at 100 M names** (the practical RAM wall for NovaSeqX-scale undetermined files).

Microbench (per record): `sha256(pickle.dumps(name))` = 0.53 µs vs builtin `hash()` = 0.022 µs; rbloom `add`/`contains` with the current sha256(pickle) hash ≈ 1.1 µs each — the hash function, not the set logic, dominates the approx path.

### Problems beyond speed

1. **The released `--method approx` CLI path is broken.** `main()` passes the lane→pair **dict** to `build_bloom_filter()`, which iterates it as if it were a list of files → `AttributeError: 'str' object has no attribute 'name'` (reproduced 2026-11-04). Approx has only ever worked via direct function calls.
2. **The two methods do opposite things.** exact purges the **undetermined** file (removes undetermined reads whose names appear in assigned files — the documented intent). approx purges the **assigned** file (removes assigned reads whose full names appear in the undetermined bloom). Output files, counts, and semantics differ.
3. **The comparison key differs per method**: exact uses `name.split(" ")[0]` (name0 = instrument:run:lane:tile:x:y), approx uses the **full** `record.name` (including the index sequence). On real Illumina data name0 is globally unique per cluster, so the difference is usually small — but on data with duplicate name0s (multi-run concatenations, some platforms) the methods diverge.
4. **Fragile concurrency**: `mp.set_start_method("fork", force=True)`, three interchangeable threading backends (`loky`/`mp_pool`/`mp_manager`), global module variables shared across forked workers, results shipped as pickled sets over process boundaries, 1 M-record intersection-update batching.
5. **Three passes over the data per undetermined lane**, and in multi-lane runs the assigned files are re-scanned **once per undetermined lane** (main loop: `for un_lane in undetermined_path: build set → scan ALL assigned → write`).
6. Bloom sized for 1 B items by default (`--max-items 1e9`): a 2.4 GB bit array even when the run has 10 M reads.

### What a rewrite buys

- 1 pass structure per concern, 1 concurrency model (threads), 0 pickling of sets.
- ~5–10× less memory (u64 hashes vs strings; ~2 GB at 100 M names vs ~10–12 GB).
- 2–5× wall time at production scale (I/O-bound floor is gzip; per-record CPU is 10–50× faster).
- Both methods work, do the same thing, with documented (or explicitly flagged) differences.
- Single static binary; no Python/deps to ship.

## 2. Goals / non-goals

**Goals**

- CLI drop-in replacement (same flags, same file-discovery rules, same output naming; see §7 checklist).
- exact + approx membership, both purging **undetermined** (canonical direction, §3).
- Correctness-first: differential-tested against the current Python exact path.
- Memory bounded and documented (§6).

**Non-goals**

- Python bindings (this is a standalone tool, not a library).
- Changing the on-disk FASTQ/BCL formats.
- Parallel decompression *within* one gzip stream (across files/lanes is enough; revisit if profiling shows a wall).
- Index-file-aware purging (R1/R2/I1/I2 collapsing stays as a file-selection rule only).

## 3. Canonical semantics (decision points)

The rewrite must decide the three divergences above. Recommendation:

| # | Decision | Recommended | Rationale |
|---|---|---|---|
| D1 | Purge direction | **Purge undetermined** in both methods (remove undetermined reads whose name appears in any assigned file) | Matches README intent and the exact method's behaviour; approx's current direction is an accident of the broken path. |
| D2 | Comparison key | **name0** (`record.name` up to first space) in both methods | name0 identifies the cluster (lane/tile/x/y); index sequences are not part of identity. exact already uses it. |
| D3 | Preserve approx's old behaviour? | Optional `--purge-assigned` flag (default off) that inverts D1, for anyone relying on it | Cheap to implement (same pipeline, swap roles); keeps the door open. |

Each decision must be called out in the changelog; D1/D2 change observable output for edge-case inputs.

## 4. Architecture

```
fastq-purge/
└── (Rust workspace)
    ├── Cargo.toml
    └── src/
        ├── main.rs        # entry point: parse → validate → dispatch
        ├── cli.rs         # clap definitions + validation (§7 parity rules)
        ├── fastq/
        │   ├── mod.rs
        │   ├── read.rs    # streaming FASTQ(.gz) reader: name0 + raw record
        │   └── write.rs   # gzip writer (level 7), single- and paired-stream
        ├── membership/
        │   ├── mod.rs
        │   ├── exact.rs   # HashSet<u64> of xxhash64(name0)
        │   └── bloom.rs   # bitvec + xxh128 (two 64-bit halves)
        └── pipeline.rs    # passes, rayon orchestration, progress reporting
```

Dependencies: `clap` (derive), `flate2` (miniz_oxide backend; `read::MultiGzDecoder` for multi-member input, `write::GzEncoder` level 7 for output), `xxhash-rust` (feature `xxh64` — the crate has **no default features**), later `rayon`, `bitvec`; `anyhow` for the bin. M1 uses a ~70-line minimal level-gated logger instead of `tracing` (added only if M3/M4 need structured logging).

Layout note (M2): the crate is a **lib** (`src/lib.rs` exposing `cli`, `discovery`, `fastq`, `logging`, `membership`, `pipeline`) plus two bins — `fastq-purge-rs` (the tool) and `fastq-dump` (dev/verification: `names`/`records` dumps for differential testing, kept for M6) — and `examples/inflate_bench.rs` (dev-only inflate throughput probe).

### Key data structures

```rust
// fastq/read.rs
pub struct NameCursor { /* position of next record in the decompressed stream */ }

/// Extract name0 (bytes up to the first b' ') from the header line,
/// then advance to the start of the next record.
pub trait NameReader {
    fn next_name0(&mut self) -> Option<Result<&[u8], FastqError>>;
}

// write.rs: paired writer that consumes (r1_bytes, r2_bytes) in lockstep.

// membership/exact.rs
pub struct ExactSet { inner: HashSet<u64> }
impl ExactSet {
    pub fn insert_name(&mut self, name0: &[u8]) {
        self.inner.insert(xxhash_h64(name0, 0));
    }
    pub fn contains_name(&self, name0: &[u8]) -> bool {
        self.inner.contains(&xxhash_h64(name0, 0))
    }
}

// membership/bloom.rs
pub struct Bloom64 {
    bits: BitVec,          // size per §6
    k: u32,                // 2 (one per 64-bit half of xxh128)
}
impl Bloom64 {
    pub fn new(n_items: u64, fpr: f64) -> Self;  // m, k per §6 formulas
    pub fn insert_name(&mut self, name0: &[u8]);  // xxh128 → two bit probes
    pub fn contains_name(&self, name0: &[u8]) -> bool;
}
```

Name extraction detail (the hot primitive): read the decompressed stream in 256 KiB chunks; find the next `\n` (header end), then the first ` ` in the header (name0 end — if absent, name0 = whole header); hash the slice; then skip 3 more `\n`s (FASTQ is strictly 4 lines/record). No per-record allocation, no full-record parse for membership passes. For the write pass, buffer the 4 raw lines as one owned `Vec<u8>` per record (only for kept records).

## 5. Pass structure

Canonical pipeline (D1/D2):

```
pass 1 — build membership          (one thread per undetermined lane, rayon)
  for each undetermined lane L:
      stream undetermined R1 (L):  name0 → insert into set/bloom for L
  (multi-lane: K independent sets, one per undetermined lane — the current
   Python re-scans assigned K times; we scan once and test against all K)

pass 2 — find duplicates           (one thread per assigned file, rayon)
  for each assigned file A (all lanes, R1-only after pair/index collapse):
      stream A: for each name0: for each undetermined lane L with a set:
          if set_L.contains(name0): dups_L.add(name0)     // or count only

pass 3 — write purged undetermined (one thread per undetermined lane, rayon)
  for each undetermined lane L:
      if dups_L empty → log "Found no duplicates..." ; write NOTHING (parity!)
      else stream (R1, R2) pair:
          if !set_of(dups_L).contains(name0): write r1, r2
```

Notes:

- Pass 2 only needs the **set of duplicate names** (not counts) — `dups_L` is a `HashSet<u64>` of hash values (same memory as pass-1 sets, typically far smaller: duplicates ≪ undetermined).
- The Python exact path ships the full duplicate set back from workers and re-tests undetermined against it; we test the same way but the "set" is in-process u64s.
- approx replaces `ExactSet` with `Bloom64` in pass 1; pass 2 membership tests become bloom queries; pass 3 tests `!bloom.contains(name0)`. False positives remove reads that should be kept (documented; `--fpr` controls it).
- Threading: `--threads N` = rayon pool size (file-level parallelism). No fork, no managers, no global state. Gzip decode cost per file is the I/O floor; 4 lanes × 2 files already fills 8 cores at typical N.
- Single-undetermined-lane + single-assigned-file (the common small case) still benefits: no process startup, no pickle, tighter inner loop.

## 6. Hash & memory math

### Hash choice: xxhash

- `xxh64` for name0 identity: ~GB/s/core, 8 B output.
- `xxh128` for the bloom: two independent 64-bit halves → k = 2 probe functions (classic optimal-k for m/n ≈ 19 bits/item is ~13; double-hashing with k=2 at 19 bits/item gives fpr ≈ (e^−19/... ) — **use the standard result**: with k=2 and m = 19.17·n bits, fpr ≈ 2^−(m/n)·... in practice compute empirically and size m with 20% margin; the test plan in §8 includes an empirical FPR check).
- Why not sha256: 30–50× slower than xxhash and provides no benefit for non-adversarial cluster names (current code uses it only because `rbloom` needs *an* int hash and `pickle` was chosen for generality — it costs 0.53 µs/record).

### 64-bit collision analysis (exact mode)

Expected number of colliding name pairs: `n(n−1) / 2^65`

| n (names) | 10 M | 100 M | 1 B |
|---|---|---|---|
| expected colliding pairs | 7×10⁻⁸ | 2.7×10⁻⁴ | 2.7×10⁻² |

At 100 M names there is a 0.027% chance of *any* spurious duplicate in a run. Acceptable for "exact" in practice; if a user needs cryptographic identity, add an optional `--strict` mode that stores name0 strings (memory per §7 table) — do not build it speculatively.

### Memory table (membership, pass-1 structure)

| Structure | Per name | 10 M | 100 M | 1 B |
|---|---|---|---|---|
| Python `set[str]` (current exact) | ~90–120 B | ~1 GB | **~10–12 GB** | infeasible |
| Python `rbloom` @ fpr 1e-4 (current approx, default 1e9 items) | 19.17 bits fixed | 24 MB | 240 MB | 2.4 GB |
| Rust `HashSet<u64>` (hashbrown, load ~0.8) | ~16 B | 0.13 GB | **~2 GB** | ~16 GB |
| Rust `Bloom64` @ fpr 1e-4 | 19.17 bits | 24 MB | **240 MB** | 2.4 GB |

Bloom sizing formula (used in `Bloom64::new`): `m = ⌈n · (−ln fpr) / (ln 2)²⌉` bits → 19.17 bits/item at fpr = 1e-4; `k = 2` (xxh128 halves) with +20% bit margin, or derive k = round(m/n · ln 2) if we switch to more probe functions.

Pass-3 `dups_L` sets: duplicates are a small fraction of undetermined in practice (that's the point of the tool); budget as a fraction of the pass-1 set, monitor in tests.

## 7. CLI parity checklist

Binary name: `fastq-purge` (or `fastq-purge-rs` during development). Flags mirror `parse_args()` (`main.py:42`) exactly:

| Flag | Type | Default | Parity requirement |
|---|---|---|---|
| `--undetermined-path` | file or dir, required | — | dir → **non-recursive** glob `*.fq*` / `*.fastq*`; files whose names contain `purged` are excluded (re-run safety) |
| `--assigned-path` | file/dir, `nargs+`, required | — | dir → **recursive** glob; dedupe; collapse paired/index files: if `re.sub(r"_[IR][0-9]_", "_", name)` was seen, skip with warning `Ignoring '<f>' as it is either a paired or index file`; lane key from `re.search(r"_L[0-9]{3}_", name)` → `L001`; files without it: warning + skip |
| `--output-path` | dir | parent of first undetermined file | create if missing (`mkdir -p` semantics) |
| `--method` | `exact` \| `approx` | `exact` | both must work (approx is currently broken) |
| `--max-items` | int | `1_000_000_000` | approx only; used for bloom sizing — **recommend keeping the default** for parity but logging the actual item count observed |
| `--fpr` | float | `1e-4` | approx only |
| `--threading-method` | `loky` \| `mp_pool` \| `mp_manager` | `mp_pool` | **accept and ignore**; one-time "accepted and ignored" warning in the pipeline phase (NOT in the validation block, to keep the validation log byte-identical to Python) |
| `--threads` | int | `1` | rayon pool size |
| `--log-level` | DEBUG…CRITICAL | `INFO` | same levels; same message content and order on **stdout** (see logging note below) |
| `--version` | | | print version |

Behavioural checklist (from `validate_args` + `main`):

- [ ] Undetermined pairing: key = 3rd `_`-separated token of the filename (e.g. `L001` in `Undetermined_S0_L001_R1_001.fastq.gz`); files sorted; a lane with a single file → `(file, None)` and only that file is purged/written.
- [ ] Output naming: `X.fastq.gz` → `X.purged.fastq.gz`; `X.fastq` → `X.purged.fastq`; written into `--output-path`.
- [ ] exact logs (same strings, same order): `Using exact method`, `Building python set from undetermined reads...`, `Number of undetermined reads: {n} ({mb:.2f} MB)`, `Found {n} duplicates in the undetermined file.`, and when zero duplicates: `Found no duplicates in the undetermined reads. Exiting...` with **no output file created**.
- [ ] approx logs: `Using bloom filter method`, `Building bloom filter...`, per-file `Processed {n} reads` every 1 M records, `Removed {n} reads`, `Kept {n} reads`, `Purged fastq file saved to '{path}'`.
- [ ] Output gzip compression level 7.
- [ ] Invalid path (neither file nor dir) → error log + exit code 1.
- [ ] `--log-level` gates DEBUG-only lines (memory usage, bloom size, item counts).
- [ ] Final `Done!` line.

**Logging contract (verified against the running Python CLI, 2026-11-05)**

The Python CLI configures `RichHandler(rich_tracebacks=True, tracebacks_show_locals=True, markup=True)` — despite `format="%(message)s"` in `basicConfig`, RichHandler renders each record itself: **stdout**, layout `HH:MM:SS LEVEL message  main.py:LINE`, level column fixed width, rich wrapping at console width, timestamp omitted when the second is unchanged, colour when TTY. The Rust CLI reproduces the *message* contract — identical messages, identical order, identical `--log-level` gating, on stdout — as plain bare lines. Deliberate deviation: no rich rendering (wall-clock timestamps are non-deterministic, `main.py:LINE` references Python source that does not exist in the rewrite, wrapping depends on console width). "Path … is neither a file nor a directory" is logged to stdout once before exit 1 (parity); errors with no Python equivalent (clean error for an empty undetermined dir, M1 pipeline bail) go to stderr.

New (documented) behaviours:

- [ ] `--purge-assigned` (D3): inverts direction; off by default.
- [ ] approx + exact produce identical outputs except for bloom false positives (test-asserted on collision-free datasets).
- [ ] Multi-undetermined-lane runs scan assigned files **once** (perf note, not an observable change except speed).

## 8. Test & rollout plan

### Milestone status

- [x] **M1 — CLI + discovery/validation (done 2026-11-05).** Rust workspace at `rust/fastq-purge/` (bin `fastq-purge-rs`, v0.3.0): clap CLI with all flags/defaults/choices; `validate_args` port (undetermined pairing, output naming/dir creation, assigned cleaning + lane grouping); fnmatch port for `pathlib.glob` single-component patterns. 17 unit tests (`cargo test`). Differential validation-log parity vs the Python CLI on the 5 M dataset: **byte-identical message blocks** for (a) single-lane gz run, (b) two-lane run exposing the stale-`lane` log quirk, (c) ungzipped `.fastq` single-read lane (`'None'` lines), (d) re-run with `.purged.*` outputs present (excluded); `--version` rc 0; usage errors rc 2 (both); invalid path → message on stdout + rc 1 (both, no duplicate print); DEBUG-only output-dir line gated identically. Deviations: empty undetermined dir → clean stderr error + rc 1 (Python: IndexError traceback); `--version` prints `fastq-purge-rs 0.3.0` (name/version differ by design); plain log rendering (see §7 note).
  - Findings: the `Assigned files:` per-file log reuses the last undetermined-lane loop variable (main.py:299) — reproduced for parity. `_R1_`/`_R2_` normalisation replaces each segment with the *literal* text `_R[12]_`, which is then re-globbed as an fnmatch class. `re.sub` consumes each matched `_[IR][0-9]_`/`_R[12]_` segment including its trailing underscore, so adjacent segments sharing an underscore are NOT both collapsed/normalised (`_R1_R2_` → `_R2_`/`R[12]_R2_`). The glob patterns `*.fq*` and `*.fastq*` are disjoint (`.fastq` does not contain the substring `.fq`); a name like `x.fq.fastq.gz` matches both and appears twice in `explode_path` output (no dedup — dedup happens downstream). Rust 1.99 removed `impl Display for PathBuf` (use `.display()`).

- [x] **M2 — FASTQ streaming + name0 hashing (done 2026-11-06).** `fastq/read.rs`: line-based port of dnaio 1.2.3's C parser over a growing 256 KiB-chunked buffer (no per-record allocation on the name0 path); `fastq/write.rs`: gzip (level 7) / plain writers, single + paired; `membership/exact.rs`: `ExactSet` = `HashSet<u64>` of xxhash64(name0); dev bin `fastq-dump` (`names`/`records`) for differential testing; crate restructured to lib + bins. 50 tests (17 M1 + 4 membership + 29 reader/writer integration, the full dnaio error/normalisation matrix). Differential on the 5 M dataset (5,000,000 records): name0 dump **byte-identical** to `record.name.split(" ")[0]` via dnaio; full record copy **byte-identical** plain and after gunzip. Single-threaded name pass: 3.6 s Rust (miniz_oxide inflate ≈ 628 MB/s, on par with C zlib ≈ 622 MB/s) vs 2.0 s Python *with pigz multi-threaded inflate* (xopen auto-detects `/usr/bin/pigz`) — the gap is thread count, to be closed in M5 (rayon).
  - Findings (dnaio 1.2.3, all verified empirically; the port reproduces them exactly):
    - A record is strictly 4 lines; header marker `@` **or** `>` (a `>`-header file is mis-parsed by dnaio itself — two-header misdetection, `qualities=None` — not reproduced; Illumina data never uses `>`).
    - CRLF is normalised to LF on read; the `+` line's content (a repeated name) is **dropped** — the writer always emits a bare `+`. `name0` = header minus the leading marker up to the first space byte (0x20); tabs are not separators; no space ⇒ whole name.
    - Error model (exact strings, 1-based physical line numbers): the header `@` check is **deferred until the sequence line has been read** (bad header + EOF ⇒ premature-EOF, not the header error); the `+` check fires when the `+` line is read; `Length of sequence and qualities differ` fires after the qual line at the qual line number; non-ASCII bytes abort as soon as the offending line is read (`Error in FASTQ file at unknown line: Non-ASCII characters found in record.`, fires before truncation errors); premature-EOF quotes the raw bytes of the incomplete record attempt in Python-`repr` style (single quotes, double if the content has a `'` and no `"`).
    - Empty file ⇒ 0 records; a final record without a trailing newline is fine; a trailing blank line is an error (it starts a new record attempt).
    - gzip: **multi-member aware** (Python's gzip module reads concatenated members; flate2's plain `GzDecoder` stops at the first member, so the reader uses `flate2::read::MultiGzDecoder` — which, like `gunzip`, errors on trailing non-gzip data). Python's gzip writer emits mtime=0, but zlib vs miniz-oxide at the same level produce different compressed bytes ⇒ `.gz` outputs are compared after decompression; plain `.fastq` outputs are byte-identical (verified: 317.7 MB vs 312.6 MB raw).
    - **Bug caught by the prefix bisect:** a buffer re-pack (memmove to front) mid-record invalidated the header/sequence line spans, so the deferred `@` check read a stale offset (one buffer position ahead) — a silent name0 corruption or a spurious `expected to start with '@'` error, manifesting only for certain file-size alignments (repro: any prefix ≥ 3000 records of the 5 M plain file). Fix: re-packing is permitted only at record boundaries (`read_line(allow_repack)`; only the header read may re-pack).
    - Parser overhead ≈ 1.1 s for 1.8 GB (≈ 1.6 GB/s): the ASCII check is a second full scan of every line and `rec` (the error-quote buffer) is appended per line; fusing the scans / making the quote lazy are the obvious M5 micro-opts if they show up in profiles.

- [x] **M3 — exact pipeline, single threaded (done 2026-11-06).** `pipeline/mod.rs` (dispatch; one-time compatibility warning when `--threading-method` ≠ `mp-pool` or `--threads` ≠ 1); `pipeline/exact.rs` — the canonical pipeline per §5: pass 1 builds one `ExactSet` per undetermined lane from the R1 name0s; pass 2 streams the assigned files **once** and tests each name0 hash against all K lane sets (Python re-scans the assigned files K times); pass 3 writes the purged undetermined files (paired = zip, stop at shorter, R1 name0 gates the pair). `--purge-assigned` (D3) inverts the direction: assigned reads whose name0 appears in the lane's undetermined set are removed, per-file outputs, warning for assigned lanes without an undetermined counterpart. `membership/exact.rs` gained `hash_name0`, `insert_hash`/`contains_hash` (one hash per name0, tested against all sets), `estimated_bytes`. 5 e2e tests (`tests/pipeline.rs`): paired lane, no-duplicates (no output files + exact log), `--purge-assigned`, single-read lane, gzipped end-to-end. 55 tests total (21 unit + 29 fastq + 5 pipeline).
  - Differential (5 M + 1.1 M mixed datasets): counts **identical** to Python (`Number of undetermined reads`, `Found {N} duplicates`, `Removed`/`Kept`); on the 1.1 M dataset (1 M real records + 100 k synthetic reads whose name0s are absent from the assigned files) the purged R1 and R2 outputs are **byte-identical after decompression** (100 k records kept each side). Wall time on the 5 M set: Rust 17.3 s vs Python 19.2 s single-threaded equivalent (Python's head start is pigz multi-threaded inflate; M5 closes the rest).
  - New log contract (documented in `pipeline/exact.rs`): the Python sequence is per-lane (build → scan *all* assigned → write, ×K); the Rust sequence builds all K sets, scans once, then writes per lane. Same messages, same counts, new order. Deliberate deviations from the Python log: (1) no `memory_usage` DEBUG lines (process/OS memory is not the same quantity and not useful for parity); (2) the `Number of undetermined reads` size figure is `set.len() × 16 B` (hashbrown entry estimate) vs Python's `sys.getsizeof` of a Python `set` — different quantity, same shape; (3) when a lane has 0 duplicates **both** `Found 0 duplicates in the undetermined file.` and `Found no duplicates in the undetermined reads. Exiting...` are printed (matching Python, which prints both too); (4) the `--threading-method`/`--threads` warning wording is new (see §7).
  - Findings (dnaio paired reader): `dnaio.open(r1, r2)` **validates that the two read names match** — `Reads are improperly paired. Read name '{n1}' in file 1 does not match '{n2}' in file 2.` — and raises `FileFormatError` at the first mismatched pair. Real Illumina paired files share read names (only the `1:N:0`/`2:N:0` suffix differs), so this never fires in production; the Rust pass-3 zip does **not** enforce it (silent zip on mismatched names). Identical behaviour on valid data; the deviation is documented, not implemented.

**Differential testing (the gate)**

1. Build a reference dataset where name0s are **globally unique** (fixes the synthetic-data caveat from profiling: the current generator's coordinate space is seed-invariant and produces birthday-collision name0s, which inflates exact-mode duplicate counts).
2. Run Python exact (threads=1) and the Rust binary on the same set; compare purged outputs record-for-record (names + full records), byte-identical when compression settings match (else record-identical).
3. Property tests: for approx, empirical FPR ≤ `--fpr` on a 10 M membership set + 10 M probes; duplicate count monotonic in `--max-items` sanity.

**Edge cases (unit/integration)**

- [ ] empty undetermined file; empty assigned dir (warning path)
- [ ] zero duplicates → no output files, exact log string
- [ ] all reads duplicated → empty (but existing) purged files
- [ ] undetermined R2 missing (`(file, None)` path)
- [ ] ungzipped `.fastq` input (output naming without `.gz`)
- [ ] multi-lane (4 lanes × R1/R2) with cross-lane duplicates
- [ ] assigned dir containing I1/I2 files (collapse + warning)
- [ ] filenames without `_L00[0-9]{2}_` (warning + skip)
- [ ] re-run with outputs present (purged files excluded from input discovery)

**Performance acceptance (5 M synthetic set, 16 cores)**

- [ ] exact end-to-end ≤ 3 s (vs 9.1 s Python)
- [ ] approx end-to-end ≤ 5 s (vs 38.4 s Python)
- [ ] 100 M-name membership ≤ 2.5 GB RSS (measured, documented)
- [ ] ≥ 2× the Python wall time at ≥ 50 M reads

**Rollout**

1. Milestones: (M1) CLI + discovery/validation port; (M2) FASTQ streaming + name0 hashing; (M3) exact pipeline single-threaded; (M4) bloom pipeline; (M5) rayon + progress; (M6) differential + perf validation. ~5–7 working days total.
2. Ship `fastq-purge-rs` side-by-side; keep the Python package for one release cycle; compare on a real run before declaring the Python version deprecated.
3. Changelog must call out D1/D2/D3 and the fixed approx path.

## 9. Risks & open questions

- **R1 — semantic divergence (D1/D2):** existing users may depend on approx's old direction/key. Mitigation: differential tests, `--purge-assigned`, explicit changelog. *Owner: decide before M1.*
- **R2 — name0 hash collisions (exact):** quantified in §6 (2.7×10⁻⁴ expected pairs at 100 M). Accept; optional `--strict` later.
- **R3 — gzip stream assumptions:** multi-member gzip (Illumina sometimes concatenates members) — `flate2` handles it; add a test fixture. CRLF line endings are not expected; fail loudly if a line lacks `\n` (corrupt input) rather than guessing.
- **R4 — output determinism:** keep writes strictly sequential per file so outputs are reproducible; only file-parallelism, never record-reordering.
- **R5 — `--max-items` default (1e9) allocates a 2.4 GB bloom for a 10 M read run.** Keep for parity, but log `bloom sized for {max_items} items, observed {n}` and consider a default change next major.
- **Open:** should pass 2 count duplicates per lane (for the `Found {n} duplicates` log) or return the set only? Counts are needed for parity logging — keep counts (cheap: one u64 per lane).