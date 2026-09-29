/// A user-supplied GenBank file. A cluster may contain more than one locus when
/// the file contains multiple GenBank records.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cluster {
    pub name: String,
    pub loci: Vec<Locus>,
}

/// One GenBank record and its coding-sequence features.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Locus {
    pub name: String,
    /// Zero-based, end-exclusive coordinate.
    pub start: usize,
    /// Zero-based, end-exclusive coordinate.
    pub end: usize,
    pub genes: Vec<Gene>,
}

/// The subset of a GenBank CDS needed by the upcoming alignment and plotting
/// stages. Raw nucleotide sequence is intentionally not kept after parsing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Gene {
    pub label: String,
    pub names: Vec<(String, String)>,
    /// Zero-based, end-exclusive coordinate covering the CDS feature.
    pub start: usize,
    pub end: usize,
    /// `1` for forward and `-1` for reverse strand.
    pub strand: i8,
    pub translation: String,
}
