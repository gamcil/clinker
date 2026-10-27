//! Shared analysis code for clinker.
//! Accepts annotation bytes rather than filesystem paths, so
//! native CLI and the browser/WASM frontend can use it.

mod align;
mod alignment;
mod analyse;
mod figure;
mod groups;
mod hierarchy;
mod layout;
mod model;
mod parse_genbank;
mod parse_gff;
mod plot_data;
mod prefilter;
mod project;
mod protein;
mod session;
mod synteny;
mod validation;

pub use align::ProteinMatch;
pub use alignment::{
    AlignmentCounts, AlignmentStats, AnalysisOptions, AnalysisOptionsError, PreparedAlignment,
    ProteinBatch, ProteinBatchAligner, ProteinBatchResult, ProteinPair, ProteinTileMatch,
};
#[cfg(not(target_arch = "wasm32"))]
pub use analyse::analyse_with_jobs;
pub use analyse::{
    AnalysisError, analyse, analyse_with_parse_options, parse_input_files,
    parse_input_files_with_options, parse_input_files_with_progress,
};
pub use figure::{ClusteringDendrogram, FigureProgress, FigureStage, PreparedFigure, RowOrder};
pub use groups::{GeneGroup, build_gene_groups, build_gene_groups_with_progress};
pub use hierarchy::ClusterOrdering;
pub use model::{
    Analysis, Cluster, Gene, GeneRef, InputFile, Link, Locus, MAX_SEQUENCE_COORDINATE,
};
pub use parse_genbank::{
    ParseError, ParseOptions, SequenceRange, parse_genbank, parse_genbank_with_options,
};
pub use parse_gff::{GffError, parse_gff3};
pub use plot_data::{
    GeneGroupDefinition, PlotCluster, PlotData, PlotGene, PlotGeneRef, PlotGroup, PlotLink,
    PlotLocus, build_plot_groups,
};
pub use project::{ClinkerProject, PROJECT_FORMAT, PROJECT_VERSION, ProjectError};
pub use session::{
    ANALYSIS_SESSION_FORMAT, ANALYSIS_SESSION_VERSION, AlignmentSettings, AnalysisSession,
    SessionError,
};
pub use synteny::ClusterPairSimilarity;
pub use validation::ValidationError;

pub(crate) use alignment::can_reach_identity;
pub use prefilter::KmerPrefilter;
pub(crate) use prefilter::cluster_kmer_candidate_pairs_with_progress;
#[cfg(test)]
use prefilter::{kmer_candidate_pairs, kmer_candidate_pairs_with_progress};
