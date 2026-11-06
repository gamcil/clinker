# Clinker benchmarks

This Snakemake workflow compares Python clinker, the published Rust CLI, and
the browser's WebAssembly alignment on the same GenBank fixtures. It reports
both end-to-end CLI time and alignment-stage time. For the alignment stage,
GenBank parsing happens before timing; candidate preparation and alignment are
timed, while figure post-processing is excluded.

The fixtures cover a small workload, many clusters, and many proteins. The
default matrix runs Python clinker `0.0.32`, clinker-cli `0.1.0` without
filtering, and clinker-cli with a `5`-residue k-mer and `8` shared-k-mer
threshold. The alignment-only comparison also includes unfiltered WASM. All
use a `0.30` identity cutoff and four jobs. Three repetitions run serially; the
summary reports median elapsed time and retained links. The k-mer run is
approximate, so link counts measure retention rather than exact parity.

## Run

Install Snakemake and Conda or Mamba. To build the WASM module, install
`wasm-pack` and the `wasm32-unknown-unknown` Rust target, then run from the
repository root:

```sh
cargo install wasm-pack
rustup target add wasm32-unknown-unknown
wasm-pack build crates/clinker-wasm --target web --out-dir ../../web/.build/pkg --out-name clinker_wasm -- --locked
```

Run the benchmark workflow:

```sh
cd benchmarks
snakemake --use-conda --conda-frontend conda --cores 4 --resources benchmark_slot=1
```

The workflow creates isolated Python, Rust, and Node environments, installs
the published packages, builds a small alignment harness against
`clinker-core`, and extracts the fixtures. The WASM stage uses Node worker
threads and the same batch alignment API as the browser; it needs no display
server and does not render a figure. Package build output and results stay
under this directory and are ignored by Git.

`results/raw.csv` contains each repetition, labeled `end_to_end` or
`alignment`. `results/summary.csv` reports per-fixture medians, speedup versus
Python, and retained links relative to unfiltered Rust. The Snakemake timing
files report full CLI command time. `results/metadata.json` records software
versions and machine details.

Edit `config.yaml` to change versions, repeats, thread count, or k-mer settings.
Change the Python package version in `python-env.yaml`. The fixture archives
and their checksums are in this directory.
