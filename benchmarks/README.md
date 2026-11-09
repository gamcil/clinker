# Clinker runtime benchmarks

Run from a repository checkout with Python 3 and Rust/Cargo installed:

```sh
python3 benchmarks/run.py
```

The first run creates a Python virtual environment, installs Python clinker,
builds the Rust CLI, and extracts the fixtures. These setup steps are outside
the timed runs. The script runs each CLI serially and writes raw repetitions,
summary CSV, metadata, logs, and link files under `benchmarks/results/`.

Run only the smallest fixture with:

```sh
python3 benchmarks/run.py --dataset bgcs_gcf_export_5
```

Options include `--dataset`, `--threads`, `--repeats`, `--identity`,
`--kmer-size`, and `--min-shared-kmers`. The default matrix uses five fixtures,
three modes (Python, Rust, Rust with the k-mer prefilter), and three
repetitions. Rust counters report potential protein pairs, pairs skipped
before alignment, and alignments performed; retained links are counted from
each output.

| Fixture | Workload | Clusters | Proteins |
| --- | --- | ---: | ---: |
| `bgcs_gcf_export_5` | Small | 15 | 221 |
| `bgcs_gcf_export_10` | Medium, balanced | 50 | 490 |
| `bgcs_gcf_export_9` | Protein-heavy | 61 | 2,110 |
| `bgcs_gcf_export_11` | Mixed | 75 | 1,359 |
| `bgcs_gcf_export_6` | Cluster-heavy | 150 | 1,424 |

Fixtures 10 and 11 are fixed samples from the existing export archives. Fixture
10 samples 50 clusters from export 6. Fixture 11 combines 50 clusters sampled
from export 6 with 25 from export 9. The archived samples and checksums make
the benchmark inputs stable across runs.

The Python environment and extracted inputs are cached in `benchmarks/.cache/`;
Cargo uses the repository's usual `target/` directory. The first run needs
network access to install Python clinker and fetch any uncached Rust crates.
On a remote node, put these caches on scratch storage if the checkout volume is
small:

```sh
CLINKER_BENCH_CACHE_DIR=/scratch/$USER/clinker-cache \
CARGO_TARGET_DIR=/scratch/$USER/clinker-target \
python3 benchmarks/run.py
```

Results, build files, and caches are ignored by Git. Fixture archives and
checksums are in `benchmarks/fixtures/` and `benchmarks/SHA256SUMS`.
