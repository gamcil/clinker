use serde::{Deserialize, Serialize};

/// Largest sequence coordinate represented exactly by the JavaScript renderer.
pub const MAX_SEQUENCE_COORDINATE: u64 = (1 << 53) - 1;

#[derive(Debug, Clone, Copy)]
pub struct InputFile<'a> {
    pub name: &'a str,
    pub bytes: &'a [u8],
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Analysis {
    pub clusters: Vec<Cluster>,
    pub links: Vec<Link>,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct Cluster {
    pub name: String,
    /// Original uploaded/input file name, retained in plot data for provenance.
    pub source_file: Option<String>,
    /// Organism parsed from the GenBank SOURCE/ORGANISM header when unambiguous.
    pub organism: Option<String>,
    pub loci: Vec<Locus>,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct Locus {
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub genes: Vec<Gene>,
}

#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct Gene {
    pub label: String,
    pub names: Vec<(String, String)>,
    pub start: u64,
    pub end: u64,
    pub strand: i8,
    pub translation: String,
}

/// A stable position within this analysis result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash, Deserialize, Serialize)]
pub struct GeneRef {
    pub cluster: usize,
    pub locus: usize,
    pub gene: usize,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Link {
    pub query: GeneRef,
    pub target: GeneRef,
    pub identity: f32,
    pub similarity: f32,
}
