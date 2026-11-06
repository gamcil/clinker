# Clinker runtime benchmarks

This Snakemake workflow compares Python clinker with the Rust CLI on shared
GenBank fixtures. It times the complete commands, including parsing, alignment,
and writing retained links. The Rust CLI also reports potential protein pairs,
pairs skipped before alignment, and alignments performed.

## Run

From a repository checkout, install Snakemake and Conda, then launch from any
directory:

```sh
./benchmarks/run.sh
```

The default run covers all three fixtures and three modes: Python, Rust without
filtering, and Rust with a 5-residue k-mer / 8 shared-kmer threshold. To run the
small fixture only:

```sh
./benchmarks/run.sh --config dataset=bgcs_gcf_export_5
```

The launcher uses four cores by default, or `SLURM_CPUS_PER_TASK` on a Slurm
allocation. Set `SNAKEMAKE_CORES` to override this. Benchmark runs are
serialized; set `--config threads=N` to match the CPUs allocated to the job.

The first run needs network access to create Conda environments, install Python
clinker, and download Rust dependencies. The Rust CLI is built from this
checkout. Tool versions and the Rust source revision are recorded in
`benchmarks/results/metadata.json`; Cargo build output is saved in
`benchmarks/results/logs/rust-install.log`.

`benchmarks/results/raw.csv` contains each repetition. `summary.csv` reports
medians, speedups, retained links, and Rust alignment-work counts. Python count
fields are blank because its CLI does not expose those counters. Results,
extracted fixtures, and environments are ignored by Git. Fixtures and their
checksums are in `benchmarks/fixtures/` and `benchmarks/SHA256SUMS`.
