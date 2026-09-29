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

/// A reusable global protein aligner.
///
/// `bio` retains its dynamic-programming and traceback buffers between calls.
/// Construct one with capacities large enough for a batch of comparisons, then
/// call [`Self::compare`] for each protein pair in that batch.
pub(crate) struct ProteinAligner {
    aligner: Aligner<fn(u8, u8) -> i32>,
}

impl ProteinAligner {
    pub(crate) fn with_capacity(query_length: usize, target_length: usize) -> Self {
        let scoring = Scoring::new(-20, -1, protein_score as fn(u8, u8) -> i32);
        Self {
            aligner: Aligner::with_capacity_and_scoring(query_length, target_length, scoring),
        }
    }

    pub(crate) fn compare(&mut self, query: &[u8], target: &[u8]) -> ProteinMatch {
        let alignment = self.aligner.global(query, target);
        metrics_from_operations(&alignment.operations, query, target)
    }
}

fn protein_score(query_amino_acid: u8, target_amino_acid: u8) -> i32 {
    2 * blosum62(query_amino_acid, target_amino_acid)
}

fn metrics_from_operations(
    operations: &[AlignmentOperation],
    query: &[u8],
    target: &[u8],
) -> ProteinMatch {
    let (mut query_index, mut target_index) = (0, 0);
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
