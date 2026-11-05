# clinker-core

Clinker compares gene clusters from GenBank or GFF3/FASTA files and prepares
editable figures. `clinker-core` provides the shared Rust parsing, alignment,
grouping, ordering, and figure-preparation library used by the native CLI and
browser app.

```rust
use clinker_core::{AnalysisOptions, InputFile, RowOrder, analyse};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bytes = std::fs::read("cluster.gbk")?;
    let inputs = [InputFile { name: "cluster.gbk", bytes: &bytes }];
    let analysis = analyse(&inputs, AnalysisOptions::default())?;
    let figure = analysis.prepare_figure(RowOrder::Similarity, &[], |_| {})?;
    println!("{} clusters", figure.plot_data.clusters.len());
    Ok(())
}
```

See the [API documentation](https://docs.rs/clinker-core) for library details,
and the [Clinker README](https://github.com/gamcil/clinker/blob/rust-port/README.md)
for the project overview, browser app, and citation.
