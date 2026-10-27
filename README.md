# clinker

Generate interactive gene cluster comparison figures from GenBank or GFF3/FASTA
files, inspect homology links, and export figures as SVG.

<img src="images/figure.png" alt="bua cluster and homologues" width="700">

## Applications

- **Web app:** load files, compare clusters, add trees and annotations, and edit
  figures locally in your browser. Input file contents are not uploaded.
- **Native CLI:** run the same Rust analysis core from the command line and
  generate plot data or a portable HTML figure. See [CLI usage](docs/cli.md).

This repository contains the Rust implementation. The original Python package and
its PyPI releases are maintained on
[`python-maintenance`](https://github.com/gamcil/clinker/tree/python-maintenance).
See [Python maintenance releases](docs/python-maintenance.md) for that process.

## Build from source

Run from the repository root:

```sh
cargo build --locked --release -p clinker-cli
./target/release/clinker examples/*.gbk --plot
./target/release/clinker examples/*.gbk --plot figure.html
./target/release/clinker examples/*.gbk --project comparison.json
```

The executable is named `clinker`. Cargo builds use the committed
viewer assets and do not require Node.js, npm, or a web build.

For the browser application, install Node.js/npm, Rust with the
`wasm32-unknown-unknown` target, and `wasm-pack`, then run:

```sh
sh web/build.sh
sh web/package.sh /tmp/clinker-web
```

Serve the packaged directory over HTTP. See [web development](web/README.md)
and [deployment](web/DEPLOYMENT.md) for details.

## Repository layout

```text
Cargo.toml, Cargo.lock   shared Rust workspace
crates/clinker-core/    parsing, alignment, grouping and plot data
crates/clinker-cli/     native command line and bundled viewer assets
crates/clinker-wasm/    browser bindings to the shared core
web/src/               browser app and workers
web/src/viewer/        shared figure viewer and project helpers
web/scripts/           build, package and asset synchronization
web/dist/              generated static application (ignored)
web/tests/             frontend checks
examples/              example inputs
bench/                 comparative analysis benchmarks
docs/                  architecture, CLI and release documentation
```

See [architecture and development](docs/architecture.md),
[project formats](docs/project-format.md), [feature status](docs/features.md),
and [release procedure](docs/releasing.md).
The native CLI embeds the shared local viewer. Bare `--plot` opens it;
`--open-project comparison.json` restores an edited project. Use `--no-browser`
to print the local URL without launching a browser. See [distribution packaging](docs/packaging.md)
for source crates, native archives and independent install checks. Portable HTML
remains a separate renderer-only export. Development and Pages deployment stay
on `rust-port` while the code base is finalized. The later [repository cutover](docs/releasing.md#repository-cutover)
covers merging to `master` and the first Rust release.

## Citation

Gilchrist, C.L.M. and Chooi, Y.-H. (2021). *clinker & clustermap.js: Automatic
generation of gene cluster comparison figures.* Bioinformatics.
[doi:10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007).
