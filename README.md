# clinker

Compare gene clusters from GenBank or GFF3/FASTA files, inspect homology links,
and create editable figures with SVG export. The web app processes files locally
in your browser; the native CLI uses the same Rust analysis core.

<img src="images/figure.png" alt="bua cluster and homologues" width="700">

## Native CLI

Build from the repository root, then open the interactive viewer:

```sh
cargo build --locked --release -p clinker-cli
./target/release/clinker examples/*.gbk --plot
```

Use `--plot figure.html` for a portable figure, or `--project comparison.json`
for an editable project. Cargo builds include the viewer and require no Node.js.
See [CLI usage](docs/cli.md) for installation and other options.

## Web app and development

See [web/README.md](web/README.md) for build, preview and test commands.
[Architecture](docs/architecture.md) describes the three Rust crates and shared
viewer. Further documentation covers [features](docs/features.md),
[project formats](docs/project-format.md), [packaging](docs/packaging.md) and
[releases](docs/releasing.md).

Rust development and Pages deployment remain on `rust-port`. The original Python
package is maintained on
[`python-maintenance`](https://github.com/gamcil/clinker/tree/python-maintenance);
see [Python releases](docs/python-maintenance.md).

## Citation

Gilchrist, C.L.M. and Chooi, Y.-H. (2021). *clinker & clustermap.js: Automatic
generation of gene cluster comparison figures.* Bioinformatics.
[doi:10.1093/bioinformatics/btab007](https://doi.org/10.1093/bioinformatics/btab007).
