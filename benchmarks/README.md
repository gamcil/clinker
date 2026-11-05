# Clinker runtime benchmarks

This Snakemake workflow compares the published Python and Rust command-line
packages on the same GenBank fixtures. It measures complete command elapsed
time, including parsing, alignment, and writing the retained-link summary.
Package installation and compilation happen before timed runs.

The three fixtures cover a small workload, many clusters, and many proteins.
The default matrix runs Python clinker `0.0.32`, clinker-cli `0.1.0` without
filtering, and clinker-cli with a `5`-residue k-mer and `8` shared-k-mer
threshold. All use a `0.30` identity cutoff and four jobs. Three repetitions
are run serially; the summary reports median elapsed time and median retained
links. The k-mer run is an approximate speed configuration, so link counts are
reported rather than treated as a parity requirement.

## Run

Install Snakemake and Conda or Mamba, then run:

```sh
cd benchmarks
snakemake --use-conda --conda-frontend conda --cores 4 --resources benchmark_slot=1
```

The first run creates isolated Python and Rust environments, installs the
published packages, and extracts the fixtures. Python and Rust installations
are isolated because both packages install a `clinker` executable. Benchmark
outputs and environments stay under this directory and are ignored by Git.

`results/raw.csv` contains each timed repetition. `results/summary.csv` reports
per-fixture medians, speedup relative to Python, and retained links relative to
the unfiltered Rust run. The Snakemake timing files record wall-clock duration
for each command. Run with `--cores 4` and the `benchmark_slot=1` resource so
timed jobs do not compete with one another.

Edit `config.yaml` to change the Rust package version, repeat count, thread
count, or k-mer settings. Change the Python version in `python-env.yaml`. The
fixtures are archived in `fixtures/`; `SHA256SUMS` records their checksums.
