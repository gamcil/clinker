#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cluster {
    pub name: String,
    /// Original uploaded/input file name, retained in plot data for provenance.
    pub source_file: Option<String>,
    /// Organism parsed from the GenBank SOURCE/ORGANISM header when unambiguous.
    pub organism: Option<String>,
    pub loci: Vec<Locus>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub genes: Vec<Gene>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub label: String,
    pub names: Vec<(String, String)>,
    pub start: usize,
    pub end: usize,
    pub strand: i8,
    pub translation: String,
}
