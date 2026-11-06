# Clinker runtime benchmarks

This Snakemake workflow compares Python clinker with the Rust CLI on the same
GenBank fixtures. It measures complete command elapsed time, including parsing,
alignment, and writing the retained-link summary. Package installation and
compilation happen before timed runs.

The three fixtures cover a small workload, many clusters, and many proteins.
The default matrix runs Python clinker `0.0.32`, Rust without filtering, and
Rust with a `5`-residue k-mer and `8` shared-kmer threshold. All use a `0.30`
identity cutoff and four jobs. Three repetitions run serially; the summary
reports median elapsed time and retained links. Rust also reports potential
cross-cluster protein pairs, pairs skipped before alignment, and alignments
performed. These counts describe work avoided; retained links remain the
result measure. The k-mer run is approximate.

## Run

Install Snakemake and Conda or Mamba, then run:

```sh
cd benchmarks
snakemake --use-conda --conda-frontend conda --cores 4 --resources benchmark_slot=1
```

The first run creates isolated Python and Rust environments, installs Python
clinker, builds the Rust CLI from this checkout, and extracts the fixtures.
Python and Rust installations are isolated because both install a `clinker`
executable. The Rust source revision is recorded in `results/metadata.json`.
Benchmark outputs and environments stay under this directory and are ignored
by Git.
If the Rust package build fails, its Cargo output is saved to
`results/logs/rust-install.log`. The Rust environment uses the current stable
compiler, whose exact version is recorded in `results/metadata.json`.

`results/raw.csv` contains each timed repetition. `results/summary.csv` reports
per-fixture medians, speedup relative to Python, retained links relative to the
unfiltered Rust run, and the Rust alignment-work counts. Python-only count
columns are blank because its CLI does not expose those counters. Snakemake
timing files record wall-clock duration for each command. Use `--cores 4` and
the `benchmark_slot=1` resource so timed jobs do not compete.

Edit `config.yaml` to change the repeat count, thread count, or k-mer settings.
Change the Python version in `python-env.yaml`. The fixtures are archived in
`fixtures/`; `SHA256SUMS` records their checksums.
