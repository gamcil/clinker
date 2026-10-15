//! Shared analysis code for clinker.
//! Accepts annotation bytes rather than filesystem paths, so
//! native CLI and the browser/WASM frontend can use it.

mod align;
mod analyse;
mod groups;
mod hierarchy;
mod layout;
mod model;
mod parse_genbank;
mod parse_gff;
mod plot_data;
mod prefilter;
mod protein;
mod session;
mod synteny;

pub use align::ProteinMatch;
#[cfg(not(target_arch = "wasm32"))]
pub use analyse::analyse_genbank_with_options_and_jobs;
pub use analyse::{
    Analysis, AnalysisError, AnalysisOptions, GeneRef, InputFile, Link, ProteinPair,
    ProteinTileMatch, analyse_cluster_pair, analyse_genbank, analyse_genbank_with_options,
    analyse_protein_pairs, analyse_protein_pairs_with_progress, analyse_protein_tile,
    analyse_protein_tile_with_progress, can_reach_identity, parse_input_files,
    parse_input_files_with_options, parse_input_files_with_options_and_progress,
    parse_input_files_with_progress,
};
pub use groups::{GeneGroup, build_gene_groups, build_gene_groups_with_progress};
pub use model::{Cluster, Gene, Locus};
pub use parse_genbank::{
    GenBankOptions, ParseError, ParseOptions, SequenceRange, parse_genbank,
    parse_genbank_with_options,
};
pub use parse_gff::{GffError, parse_gff3};
pub use plot_data::{GeneGroupDefinition, PlotCluster, PlotData, PlotGroup, build_plot_groups};
pub use session::{
    ANALYSIS_SESSION_FORMAT, ANALYSIS_SESSION_VERSION, AnalysisSession, SessionError,
};
pub use synteny::ClusterPairSimilarity;

pub use prefilter::{
    KmerPrefilter, cluster_kmer_candidate_pairs_with_progress, kmer_candidate_pairs,
    kmer_candidate_pairs_with_progress,
};
