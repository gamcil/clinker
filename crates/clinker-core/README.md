# clinker-core

Shared Rust parsing, protein alignment, grouping, cluster ordering and figure
preparation for the clinker CLI and browser app. Inputs are named byte slices;
no filesystem or frontend tools are required.

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

`prepare_figure` scores, orders and arranges a figure. `to_plot_data(order,
definitions)` converts existing analysis without scoring or changing locus
placement. `AnalysisSession` retains analysis for continuation; `ClinkerProject`
also preserves the editable figure.

See the [repository](https://github.com/gamcil/clinker) for usage and documentation.
