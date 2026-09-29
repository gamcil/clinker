#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Cluster {
    pub name: String,
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
