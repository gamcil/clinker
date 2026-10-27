# clinker-cli

Native clinker analysis and an embedded interactive gene cluster viewer.
The executable is named `clinker`.

After the Rust crates are published:

```sh
cargo install --locked clinker-cli
```

From a source checkout:

```sh
cargo install --locked --path crates/clinker-cli
```

```sh
clinker clusters/*.gbk --plot
clinker cluster.gff3 --plot figure.html
clinker clusters/*.gbk --project comparison.json
clinker --open-project comparison.json
```

Bare `--plot` opens the local viewer and keeps the CLI running until Ctrl+C.
`--no-browser` prints the loopback URL without launching a browser. A filename
after `--plot` writes renderer JSON or a portable, renderer-only HTML figure.
Projects retain analysis data, annotations, trees, appearance and figure edits.
GFF3 inputs need embedded FASTA or a same-stem companion FASTA file.

The complete viewer and licence notices ship inside this source package and are
embedded by a Rust-only build script. Cargo builds/installations never invoke
Node.js, npm or wasm-pack, and never download frontend assets. The viewer renders
native results without WASM analysis. Third-party viewer licences are included
under `assets/`.

See the [clinker repository](https://github.com/gamcil/clinker) for release
archives, detailed usage, development documentation and citation information.
