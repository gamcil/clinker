# Clinker runtime benchmarks

Run from a repository checkout with Python 3 and Rust/Cargo installed:

```sh
python3 benchmarks/run.py
```

The first run creates a Python virtual environment, installs Python clinker,
builds the Rust CLI, and extracts the fixtures. These setup steps are outside
the timed runs. The script runs each CLI serially and writes raw repetitions,
summary CSV, metadata, logs, and link files under `benchmarks/results/`.

Run only the small fixture with:

```sh
python3 benchmarks/run.py --dataset bgcs_gcf_export_5
```

Options include `--threads`, `--repeats`, `--identity`, `--kmer-size`, and
`--min-shared-kmers`. The default matrix uses three fixtures, three modes
(Python, Rust, Rust with the k-mer prefilter), and three repetitions. Rust
counters report potential protein pairs, pairs skipped before alignment, and
alignments performed; retained links are counted from each output.

The Python environment and extracted inputs are cached in `benchmarks/.cache/`;
Cargo uses the repository's usual `target/` directory. The first run needs
network access to install Python clinker and fetch any uncached Rust crates.
Results, build files, and caches are ignored by Git. Fixture archives and
checksums are in `benchmarks/fixtures/` and `benchmarks/SHA256SUMS`.
