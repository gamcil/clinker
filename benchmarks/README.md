# Clinker runtime benchmarks

This workflow compares Python Clinker 0.0.32, Rust clinker-cli 0.1.0, and the
browser app built from the image's source checkout. Each fixture is run with
Python, Rust, and WebAssembly; Rust and WebAssembly are also run with the k-mer
prefilter (size 5, minimum 8 shared k-mers). Runs are sequential, with up to
four threads each. The five fixtures and their checksums are in `fixtures.tsv`.

## Run

After the image has been published, download it with Apptainer and run it on a
node with a writable working directory:

```sh
apptainer pull clinker-benchmarks.sif \
  docker://ghcr.io/gamcil/clinker-benchmarks:0.1.0
mkdir -p benchmark-run
apptainer run --bind "$PWD/benchmark-run:/work" \
  clinker-benchmarks.sif --threads 4 --repeats 3
```

To publish an image from `rust-port`, tag the desired commit and push the tag:

```sh
git tag -a benchmark-v0.1.0 -m "Benchmark image 0.1.0"
git push origin benchmark-v0.1.0
```

The benchmark image workflow builds `linux/amd64` and publishes tag `0.1.0` to
GHCR. Its Actions summary records the source commit and image digest. The
workflow also supports manual dispatch once it exists on the default branch.

To build the image locally, from the repository root:

```sh
docker build -f benchmarks/Dockerfile -t clinker-benchmarks:local .
mkdir -p benchmark-run
docker run --rm -v "$PWD/benchmark-run:/work" \
  clinker-benchmarks:local --threads 4 --repeats 3
```

Use `--fixture bgcs_gcf_export_5` for the
smallest fixture; repeat `--fixture` to select several. Other options are
`--identity`, `--kmer-size`, and `--min-shared-kmers`. Run the same command
again to resume; keep `benchmark-run/` so Snakemake retains its metadata and
completed outputs. The packaged build fingerprint and settings are tracked so
changing them invalidates affected measurements. For a release
run, set `BENCHMARK_IMAGE` through Apptainer's `--env` option to the immutable
image digest so it is recorded with the results.

Results are written to `benchmark-run/results/`: `runs.csv` contains individual
measurements, `summary.csv` contains medians, and `run-metadata.json` records
software and machine details. Open `benchmark-results.html` for the table and
plots. Prepared figure JSON and logs are retained for each run.

Each measurement starts before input files are read and ends after figure data
is serialized to JSON. It includes parsing, alignment, gene grouping and
figure preparation. Process and browser startup, loading Python and JavaScript
modules, browser rendering, and writing the resulting JSON are outside the
timed interval. The browser measurement still includes its worker and WASM
initialization and transfers. Implementations use their own figure preparation
and JSON representations, which are not identical; this is a comparison of the
same workflow boundary rather than identical instructions or output bytes.

The native adapter uses the published `clinker-core` 0.1.0; the browser uses
the source checkout recorded in the image. Python uses its plotting data API,
not HTML generation. Runs report retained link and group counts. After timing,
shared links and identical gene groups are compared with unfiltered Rust using
cluster, locus and gene labels. `runs.csv` flags duplicate labels or differing
gene identities, which can make those comparisons incomplete. Link percentages
in the report compare counts, not link identities. Use a fresh working directory
when changing benchmark versions so older measurements cannot be mixed in.
