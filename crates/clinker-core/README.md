# clinker-core

Shared Rust analysis for clinker gene cluster comparison figures. The core parses
GenBank and GFF3/FASTA bytes, compares proteins, groups homologues, orders clusters
and prepares renderer data. Native and WebAssembly frontends use this same crate.

```rust
use clinker_core::{AnalysisOptions, InputFile, RowOrder, analyse_genbank};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let bytes = std::fs::read("cluster.gbk")?;
    let inputs = [InputFile { name: "cluster.gbk", bytes: &bytes }];
    let analysis = analyse_genbank(&inputs, AnalysisOptions::default())?;
    let figure = analysis.prepare_figure(RowOrder::Similarity, &[], |_| {})?;
    println!("{} clusters", figure.plot_data.clusters.len());
    Ok(())
}
```

`AnalysisSession` retains analysis state for continuation. `ClinkerProject`
combines that state with an editable figure; renderer-specific data remains JSON
so it can preserve saved annotations, trees and editor state.

The package includes its own regression fixtures. Building it needs Rust and the
target platform's native toolchain; it does not require Python, Node.js, npm,
WASM tooling or frontend assets.

See the [clinker repository](https://github.com/gamcil/clinker) for applications,
development documentation, examples and citation information.
