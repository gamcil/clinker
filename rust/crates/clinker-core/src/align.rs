use block_aligner::{
    cigar::{Cigar, Operation},
    scan_block::{Block, PaddedBytes},
    scores::{AAMatrix, BLOSUM62, Gaps, Matrix},
};

/// The similarity metrics clinker stores for one pair of proteins.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinMatch {
    pub identity: f32,
    pub similarity: f32,
}

const MIN_BLOCK_SIZE: usize = 32;
const MAX_BLOCK_SIZE: usize = 256;
const GAPS: Gaps = Gaps {
    // Doubled from clinker’s -10 / -0.5 configuration so the half-point
    // extension penalty is representable with block-aligner's integer scores.
    open: -20,
    extend: -1,
};

/// A reusable SIMD global protein aligner.
///
/// `Block`, padded sequences, and traceback storage are allocated once for a
/// tile and reused for every comparison in that tile. The adaptive 32–256
/// block range is block-aligner’s recommended setting for proteins.
pub(crate) struct ProteinAligner {
    matrix: AAMatrix,
    block: Block<true, false>,
    query: PaddedBytes,
    target: PaddedBytes,
    cigar: Cigar,
    query_bytes: Vec<u8>,
    target_bytes: Vec<u8>,
}

impl ProteinAligner {
    pub(crate) fn with_capacity(query_length: usize, target_length: usize) -> Self {
        Self {
            matrix: doubled_blosum62(),
            block: Block::new(query_length, target_length, MAX_BLOCK_SIZE),
            query: PaddedBytes::new::<AAMatrix>(query_length, MAX_BLOCK_SIZE),
            target: PaddedBytes::new::<AAMatrix>(target_length, MAX_BLOCK_SIZE),
            cigar: Cigar::new(query_length, target_length),
            query_bytes: Vec::with_capacity(query_length),
            target_bytes: Vec::with_capacity(target_length),
        }
    }

    pub(crate) fn compare(&mut self, query: &[u8], target: &[u8]) -> ProteinMatch {
        normalise_protein(query, &mut self.query_bytes);
        normalise_protein(target, &mut self.target_bytes);
        self.query
            .set_bytes::<AAMatrix>(&self.query_bytes, MAX_BLOCK_SIZE);
        self.target
            .set_bytes::<AAMatrix>(&self.target_bytes, MAX_BLOCK_SIZE);
        self.block.align(
            &self.query,
            &self.target,
            &self.matrix,
            GAPS,
            MIN_BLOCK_SIZE..=MAX_BLOCK_SIZE,
            0,
        );
        let result = self.block.res();
        self.block.trace().cigar_eq(
            &self.query,
            &self.target,
            result.query_idx,
            result.reference_idx,
            &mut self.cigar,
        );
        metrics_from_cigar(&self.cigar, query, target)
    }
}

fn doubled_blosum62() -> AAMatrix {
    let mut matrix = AAMatrix::new();
    for one in b'A'..=b'Z' {
        for two in one..=b'Z' {
            let score = BLOSUM62.get(one, two);
            if score != i8::MIN {
                matrix.set(one, two, 2 * score);
            }
        }
    }
    matrix
}

fn normalise_protein(protein: &[u8], output: &mut Vec<u8>) {
    output.clear();
    output.extend(protein.iter().map(|&amino_acid| {
        let amino_acid = amino_acid.to_ascii_uppercase();
        if b"ARNDCQEGHILKMFPSTWYVBZX".contains(&amino_acid) {
            amino_acid
        } else {
            b'X'
        }
    }));
}

fn metrics_from_cigar(cigar: &Cigar, query: &[u8], target: &[u8]) -> ProteinMatch {
    let (mut query_index, mut target_index) = (0, 0);
    let (mut alignment_columns, mut identical, mut similar) = (0_u32, 0_u32, 0_u32);

    for index in 0..cigar.len() {
        let operation = cigar.get(index);
        alignment_columns += operation.len as u32;
        match operation.op {
            Operation::M | Operation::Eq | Operation::X => {
                for _ in 0..operation.len {
                    if query[query_index] == target[target_index] {
                        identical += 1;
                    } else {
                        similar += amino_acids_are_similar(query[query_index], target[target_index])
                            as u32;
                    }
                    query_index += 1;
                    target_index += 1;
                }
            }
            Operation::I => query_index += operation.len,
            Operation::D => target_index += operation.len,
            Operation::Sentinel => unreachable!("traceback CIGAR has no sentinel operation"),
        }
    }

    ProteinMatch {
        identity: identical as f32 / alignment_columns as f32,
        similarity: (identical + similar) as f32 / alignment_columns as f32,
    }
}

fn amino_acids_are_similar(one: u8, two: u8) -> bool {
    const GROUPS: [&[u8]; 7] = [b"GAVLI", b"FYW", b"CM", b"ST", b"KRH", b"DENQ", b"P"];
    GROUPS
        .iter()
        .any(|group| group.contains(&one) && group.contains(&two))
}

#[cfg(test)]
mod tests {
    use super::ProteinAligner;

    #[test]
    fn calculates_identity_and_similarity_from_global_alignment() {
        let mut aligner = ProteinAligner::with_capacity(4, 4);

        let identical = aligner.compare(b"MST", b"MST");
        assert_eq!(identical.identity, 1.0);
        assert_eq!(identical.similarity, 1.0);

        let conservative_substitution = aligner.compare(b"MST", b"MTT");
        assert!((conservative_substitution.identity - 2.0 / 3.0).abs() < f32::EPSILON);
        assert_eq!(conservative_substitution.similarity, 1.0);

        let terminal_gap = aligner.compare(b"MST", b"MSTT");
        assert_eq!(terminal_gap.identity, 0.75);
        assert_eq!(terminal_gap.similarity, 0.75);
    }
}
