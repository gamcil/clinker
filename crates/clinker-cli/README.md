# clinker

Compare gene clusters from GenBank or GFF3/FASTA files, inspect homology links,
and create editable figures. `clinker-cli` is the native command-line
application; it uses the same Rust analysis core as the browser app.

Install the command-line application:

```sh
cargo install --locked clinker-cli
```

Open an interactive viewer or save a portable figure or editable project:

```sh
clinker clusters/*.gbk --plot
clinker cluster.gff3 --plot figure.html
clinker clusters/*.gbk --project comparison.json
```

See the [CLI guide](https://github.com/gamcil/clinker/blob/rust-port/docs/cli.md)
for options and supported inputs. The [Clinker README](https://github.com/gamcil/clinker/blob/rust-port/README.md)
introduces the project, browser app, and citation.
