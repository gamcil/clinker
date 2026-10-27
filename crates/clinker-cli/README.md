# clinker-cli

Native gene-cluster comparison and an embedded interactive viewer. Install from
a source checkout:

```sh
cargo install --locked --path crates/clinker-cli
clinker clusters/*.gbk --plot
clinker cluster.gff3 --plot figure.html
clinker clusters/*.gbk --project comparison.json
clinker --open-project comparison.json
```

Bare `--plot` opens the viewer until Ctrl+C; `--no-browser` prints its local URL.
A `.html` or `.json` filename writes a portable figure or raw renderer data.
Projects preserve analysis, annotations, trees and figure edits. GFF3 inputs
require embedded FASTA or a same-stem companion file.

Cargo builds embed the included viewer without Node.js or frontend downloads.
Third-party licences are under `assets/`. After registry publication, install
with `cargo install --locked clinker-cli`.

See the [repository CLI guide](https://github.com/gamcil/clinker/blob/rust-port/docs/cli.md)
for all options and supported platforms.
