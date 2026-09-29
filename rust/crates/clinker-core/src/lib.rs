//! Shared analysis code for clinker.
//! Accepts GenBank bytes rather than filesystem paths, so
//! native CLI and the browser/WASM frontend can use it.

mod align;
mod analyse;
mod model;
mod parse_genbank;
mod plot_data;

pub use align::ProteinMatch;
pub use analyse::{
    Analysis, AnalysisError, AnalysisOptions, GeneRef, InputFile, Link, ProteinTileMatch,
    analyse_cluster_pair, analyse_genbank, analyse_protein_tile, parse_input_files,
};
pub use model::{Cluster, Gene, Locus};
pub use parse_genbank::{ParseError, parse_genbank};
pub use plot_data::PlotData;
