# Alignment runtime benchmark

Compare the legacy Python implementation in `clinker/`, optimized native Rust,
and the real browser WASM/Web Worker implementation. The benchmark stops after
retained homology links are collected: it does not run grouping, ordering,
layout, rendering, or save large link tables.

Default settings: both local datasets (`bgcs/`, `bgcs_gcf_export_1/`), five
iterations, four alignment workers, 30% identity, seed 1. Each iteration contains
all ten dataset/configuration combinations in a seeded shuffled order: 50 runs.

| Configuration | Implementation | K-mer filter |
| --- | --- | --- |
| python | Legacy Python source in this checkout | Off |
| native | Release-mode clinker-core | Off |
| wasm | Chromium Web Workers and production WASM | Off |
| native-prefilter | Release-mode clinker-core | k=3, minimum count=3 |
| wasm-prefilter | Chromium Web Workers and production WASM | k=3, minimum count=3 |

Both native configurations use the CLI's shared core workflow, via
`Analysis::append_clusters_with_options_and_jobs` on an empty analysis.
Native and browser runs share the core batch planner and comparison loop.
Batches contain at most 1,600 pairs; Rayon and Web Workers handle dispatch on
their respective platforms. Prefiltering changes candidate selection, and both
modes reuse alignment buffers across batches. The CLI exposes the same settings with
`--prefilter --kmer-size K --min-shared-kmers COUNT`.
The browser imports `rust/web/alignment.mjs`, shared with the app, and uses its
production tiling, workers and result collection with an explicit worker limit.
No implementation substitutes a new alignment algorithm.

## Build

From the repository root, use Python 3.9+ with the legacy dependencies. The
existing local `vnv/bin/python` works on this development machine. To create a
separate environment elsewhere:

```sh
python3 -m venv bench/.venv
bench/.venv/bin/python -m pip install -e .
```

The adapters import the Python source from this checkout, not another installed
clinker package. Package versions, Python source hashes, git status, native
binary hash and WASM hash are recorded. For publication, use committed source
revisions and retain the recorded Python package versions.

```sh
cargo build --manifest-path rust/Cargo.toml --locked --release -p clinker-core --example benchmark
sh rust/web/build.sh
cd rust/web && npx playwright install chromium && cd ../..
```

Node and the web build's WASM prerequisites are required. Run on macOS/Linux
(the timeout cleanup uses process groups), with no other benchmarking jobs
running. Chromium is headless: this tests CPU/WASM analysis, not rendering.
Browser launch and WASM compilation are outside the processing timer.

## Provided datasets and the strict input check

Before timing, all three implementations parse each entire selected dataset.
The runner requires identical protein sequences, including duplicates and
per-file membership. Within-file ordering is ignored because circular-origin
rotation can reorder otherwise identical genes. A mismatch aborts all timing
and saves the parsed sequences for diagnosis; there is no silently permissive
mode.

The original local datasets have three CDS features whose fuzzy coordinates
fall outside a fetched region at its boundary. Python retains their supplied
translations, while Rust excludes them:

- `bgcs/KB908593.1_492644-545253.gbk`: `EOA86445.1`
- `bgcs/MZNU01000266.1_149617-201695.gbk`: `OWP01624.1`
- `bgcs_gcf_export_1/PDHQ01000005.1_3063-46772.gbk`: `RDA95433.1`

For these datasets, use **`--prepare-inputs`**. It makes separate copies under
the output directory. A fixed coordinate-based rule excludes a CDS when any
part has a fuzzy end before the first base, or a fuzzy start after the last
base, represented by the one-base boundary fragments in these downloaded files.
The rule does not inspect timing or alignment results. It records every
removal in `preparation.json`; originals and ordinary partial CDS features
elsewhere are left alone. Unchanged files are copied byte-for-byte; changed
records are written by Biopython. All engines receive those exact same copies.
The original and prepared file checksums are both retained.

After preparation, the complete datasets contain 204 clusters/2,019 proteins
and 110 clusters/1,594 proteins, respectively. Exact three-way sequence equality
has been verified. Without preparation, raw-input mode correctly stops on the
mismatches. These are benchmarks of the explicitly prepared datasets, not an
unqualified comparison of parsing all original annotations identically.

## Run

Use the dependency-equipped interpreter (`vnv/bin/python` locally, or replace it
with `bench/.venv/bin/python`).

```sh
# Full requested experiment: 5 iterations × 2 datasets × 5 configurations.
vnv/bin/python bench/run.py --prepare-inputs --output benchmark-results/full

# Check the complete input sets without running any alignments.
vnv/bin/python bench/run.py --prepare-inputs --validate-only --output benchmark-results/validation

# Fast integration check: two smallest files per dataset, one iteration.
# This deliberately uses a small subset, not the publication workload.
vnv/bin/python bench/run.py --smoke --output benchmark-results/smoke

# Rust only, one iteration on each full dataset, with a 5/8 prefilter.
vnv/bin/python bench/run.py --prepare-inputs --repeats 1 \
  --configurations native wasm native-prefilter wasm-prefilter \
  --kmer-size 5 --min-shared-kmers 8 --output benchmark-results/rust-5-8
```

Output directories must be empty. Omitting `--output` creates a timestamped
directory. Options include `--datasets DIR [DIR ...]`, `--repeats 5`,
`--workers 4`, `--identity 0.3`, `--seed 1`, `--timeout 3600`, and
`--python /path/to/python`. Timeout is per adapter invocation, in seconds.
`--configurations` selects the measured implementations; protein preflight
still checks all three. `--kmer-size` and `--min-shared-kmers` control both Rust
prefilter configurations and are recorded in the run CSV and metadata.
Large Python runs can take a long time; increase it if needed. Failure/timeout
saves the partial CSV and logs, terminates the run and does not silently retry,
discard, or report a partial run as complete. Run one experiment at a time.

## What is timed

All engines stage input bytes before starting their monotonic processing clock.
The Python adapter uses the same Biopython `SeqIO`/`Locus.from_seqrecord`/`Cluster`
construction as `cluster_from_genbank`, reading staged bytes in memory rather
than reopening paths. Native and WASM call their existing parsers. Decoding,
parsing and translation construction are inside `parse_ms`.

`total_ms` is the primary measurement: elapsed time from parsing start to
collected homology links, including optional prefiltering, alignment-pool
construction and coordination. Workload counts come from shared core preparation
and execution, without repeating candidate selection. Reporting happens after timing.

Some adapters also emit phase timings as diagnostics. Their boundaries differ
between implementations and should not be compared. Native reports only parse
and total time; filtering and alignment use one shared core call.

Module imports, process/browser startup, file staging, WASM compilation,
validation, post-run diagnostics, final JSON/CSV writing and teardown are excluded.
The parser worker is initialized before timing; the alignment pool is newly
created inside the alignment phase for every measured run in all engines.

Use `total_ms` as the primary result and phase columns as explanation. This is
a staged-input processing benchmark, not command-launch-to-exit or a browser
upload/UI benchmark. Worker count is configured parallelism, not identical total
OS threads; browser worker count may be lower if fewer tiles exist. Native and
Python use their existing task granularity, which is intentionally different.
BLAS/OpenMP helper threads are limited to one in subprocesses.

## Work counts and output

The k-mer-off Rust paths still perform their production length-ratio rejection;
Python v1 does not. Thus matched inputs/settings do not imply identical numbers
of aligner calls. `possible_pairs`, `candidate_pairs`, `aligned_pairs` and
`retained_links` make this explicit. Workload counts are collected by the core workflow and reported after timing;
`eligible_pairs` is unavailable in the current adapters. The runner verifies matching
native/WASM unfiltered execution counts and filtered candidate counts across
each iteration; alignment results are not required to match Python numerically.

- `runs.csv`: each measured run, all phase timings, settings, counts and status.
- `summary.csv`: per-dataset/configuration median, minimum and maximum total time.
- `datasets.json`: actual file order and checksums; inputs checked again at end.
- `source-datasets.json`, `preparation.json`, `inputs/`: original provenance and
  prepared copies when requested.
- `validation.json`, `proteins-*.json`: shared input counts and exact sequences.
- `schedule.json`, `metadata.json`, stderr logs: execution order and provenance.

No filesystem-cache flushing is attempted. Preflight and staging warm the input
cache; processing uses fresh analysis state/pools each time. Use the same machine
and power settings, close unrelated jobs, keep all runs, and report variation
across repetitions. The harness cannot eliminate all OS background activity.

For implementation verification, run the smoke command, `npm --prefix rust/web
test`, the normal browser workflows, and `cargo clippy --manifest-path
rust/Cargo.toml --locked -p clinker-core --example benchmark -- -D warnings`.
