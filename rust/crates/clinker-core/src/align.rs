use bio::{
    alignment::{
        AlignmentOperation,
        pairwise::{Aligner, Scoring},
    },
    scores::blosum62,
};

/// The similarity metrics clinker stores for one pair of proteins.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinMatch {
    pub identity: f32,
    pub similarity: f32,
}

/// Globally align two protein sequences with scaled BLOSUM62 affine scoring.
///
/// Scores are multiplied by two so integer scoring preserves the original
/// Python defaults of gap open -10 and gap extension -0.5.
pub fn compare_proteins(query: &[u8], target: &[u8]) -> ProteinMatch {
    let scoring = Scoring::new(-20, -1, |query_amino_acid, target_amino_acid| {
        2 * blosum62(query_amino_acid, target_amino_acid)
    });
    let mut aligner = Aligner::with_capacity_and_scoring(query.len(), target.len(), scoring);
    let alignment = aligner.global(query, target);
    metrics_from_operations(&alignment.operations, query, target)
}

fn metrics_from_operations(
    operations: &[AlignmentOperation],
    query: &[u8],
    target: &[u8],
) -> ProteinMatch {
    let (mut query_index, mut target_index) = (0, 0);
    // This deliberately counts gap columns. The original clinker
    // `compute_identity` divides by the rendered length of a global
    // Biopython alignment, rather than its number of residue-residue columns.
    // Without this, a short matching segment surrounded by long gaps can pass
    // the homology cutoff and collapse unrelated genes into one group.
    let (mut alignment_columns, mut identical, mut similar) = (0_u32, 0_u32, 0_u32);
    for operation in operations {
        alignment_columns += 1;
        match operation {
            AlignmentOperation::Match => {
                identical += 1;
                query_index += 1;
                target_index += 1;
            }
            AlignmentOperation::Subst => {
                similar += amino_acids_are_similar(query[query_index], target[target_index]) as u32;
                query_index += 1;
                target_index += 1;
            }
            // rust-bio's names are relative to the first (query/x) sequence:
            // an insertion consumes x, while a deletion consumes y.
            AlignmentOperation::Ins => query_index += 1,
            AlignmentOperation::Del => target_index += 1,
            _ => unreachable!("global alignment does not contain clipping operations"),
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
    use super::compare_proteins;

    #[test]
    fn calculates_identity_and_similarity_from_global_alignment() {
        let identical = compare_proteins(b"MST", b"MST");
        assert_eq!(identical.identity, 1.0);
        assert_eq!(identical.similarity, 1.0);

        let conservative_substitution = compare_proteins(b"MST", b"MTT");
        assert!((conservative_substitution.identity - 2.0 / 3.0).abs() < f32::EPSILON);
        assert_eq!(conservative_substitution.similarity, 1.0);

        let terminal_gap = compare_proteins(b"MST", b"MSTT");
        assert_eq!(terminal_gap.identity, 0.75);
        assert_eq!(terminal_gap.similarity, 0.75);
    }
}
