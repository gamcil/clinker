//! Shared analysis code for clinker.
//! Accepts GenBank bytes rather than filesystem paths, so
//! native CLI and the browser/WASM frontend can use it.

mod align;
mod analyse;
mod groups;
mod hierarchy;
mod layout;
mod model;
mod parse_genbank;
mod plot_data;
mod synteny;

pub use align::ProteinMatch;
pub use analyse::{
    Analysis, AnalysisError, AnalysisOptions, GeneRef, InputFile, KmerPrefilter, Link, ProteinPair,
    ProteinTileMatch, analyse_cluster_pair, analyse_genbank, analyse_protein_pairs,
    analyse_protein_pairs_with_progress, analyse_protein_tile, analyse_protein_tile_with_progress,
    can_reach_identity, kmer_candidate_pairs, parse_input_files,
};
pub use groups::{GeneGroup, build_gene_groups};
pub use model::{Cluster, Gene, Locus};
pub use parse_genbank::{ParseError, parse_genbank};
pub use plot_data::{PlotCluster, PlotData, PlotGroup};
pub use synteny::DEFAULT_CONTIGUITY_WEIGHT;
