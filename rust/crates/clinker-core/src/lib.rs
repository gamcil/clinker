//! Shared, platform-independent analysis code for clinker.
//!
//! This crate deliberately accepts GenBank bytes rather than filesystem paths,
//! so both the native CLI and the browser/WASM frontend can use it.

mod align;
mod analyse;
mod model;
mod parse_genbank;
mod plot_data;

pub use align::{ProteinMatch, compare_proteins};
pub use analyse::{
    Analysis, AnalysisError, AnalysisOptions, GeneRef, InputFile, Link, analyse_genbank,
};
pub use model::{Cluster, Gene, Locus};
pub use parse_genbank::{ParseError, parse_genbank};
pub use plot_data::PlotData;
