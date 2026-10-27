# clinker-core

Shared Rust analysis for clinker gene cluster comparison figures. The core parses
GenBank and GFF3/FASTA bytes, compares proteins, groups homologues, orders clusters
and prepares renderer data. Native and WebAssembly frontends use this same crate.

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

`analyse` accepts GenBank and GFF3/FASTA inputs. Use
`analyse_with_parse_options` for explicit parsing settings, or the native
`analyse_with_jobs` to select concurrency. Parsing APIs share `ParseOptions`;
`parse_input_files_with_progress` takes those options and a progress callback.

Analysis entry points and worker preparation reject invalid identity cutoffs or
prefilter settings with `AnalysisOptionsError`. Serial and native parallel
append methods both take `AnalysisOptions` and return work statistics. Worker
consumers use `PreparedAlignment` and `ProteinBatchAligner`; candidate-selection
helpers remain internal.

Scoring and ordering return validation errors for malformed analysis data.
Append operations validate existing and incoming data before changing the
analysis. Caller-supplied matrices are kept out of the public ordering API;
`prepare_figure` returns the computed scores, ordering and hierarchy together.

`AnalysisSession` retains analysis state for continuation. `ClinkerProject`
combines that state with an editable figure; renderer-specific data remains JSON
so it can preserve saved annotations, trees and editor state.

The package includes its own regression fixtures. Building it needs Rust and the
target platform's native toolchain; it does not require Python, Node.js, npm,
WASM tooling or frontend assets.

See the [clinker repository](https://github.com/gamcil/clinker) for applications,
development documentation, examples and citation information.
