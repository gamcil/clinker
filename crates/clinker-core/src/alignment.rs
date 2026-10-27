use std::ops::Range;

use serde::{Deserialize, Serialize};
use thiserror::Error;

use crate::{
    Cluster, GeneRef, KmerPrefilter, Link, align::ProteinAligner,
    cluster_kmer_candidate_pairs_with_progress,
};

#[derive(Debug, Clone, Copy)]
pub struct AnalysisOptions {
    pub identity_cutoff: f32,
    /// Optional approximate candidate selection before alignment.
    pub prefilter: Option<KmerPrefilter>,
}

#[derive(Debug, Error)]
pub enum AnalysisOptionsError {
    #[error("identity cutoff must be between 0 and 1, got {0}")]
    IdentityCutoff(f32),
    #[error("prefilter requires k-mer size 1..=8 and a positive shared-word count")]
    Prefilter,
}

impl AnalysisOptions {
    pub fn validate(self) -> Result<(), AnalysisOptionsError> {
        if !self.identity_cutoff.is_finite() || !(0.0..=1.0).contains(&self.identity_cutoff) {
            return Err(AnalysisOptionsError::IdentityCutoff(self.identity_cutoff));
        }
        if let Some(filter) = self.prefilter
            && (!(1..=8).contains(&filter.kmer_size) || filter.min_shared_kmers == 0)
        {
            return Err(AnalysisOptionsError::Prefilter);
        }
        Ok(())
    }
}

/// Work performed by the alignment workflow, after length and optional k-mer filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentStats {
    pub aligned_pairs: usize,
    pub counts: AlignmentCounts,
}

impl Default for AnalysisOptions {
    fn default() -> Self {
        Self {
            identity_cutoff: 0.30,
            prefilter: None,
        }
    }
}

/// A retained match between two positions in a protein-comparison tile.
///
/// The indices refer to the caller-provided query and target slices rather
/// than loci or clusters. Browser workers use this compact result while the
/// coordinator retains the corresponding gene references.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinTileMatch {
    pub query_index: usize,
    pub target_index: usize,
    pub identity: f32,
    pub similarity: f32,
}

/// A pair of positions in caller-provided protein slices.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ProteinPair {
    pub query_index: usize,
    pub target_index: usize,
}

/// The largest possible global-alignment identity is the shorter sequence
/// divided by the longer one. Rejecting pairs below this bound is exact.
pub(crate) fn can_reach_identity(
    query_length: usize,
    target_length: usize,
    identity_cutoff: f32,
) -> bool {
    if identity_cutoff <= 0.0 {
        return true;
    }
    let longer = query_length.max(target_length);
    longer != 0 && query_length.min(target_length) as f32 / longer as f32 >= identity_cutoff
}

const TILE_SIDE: usize = 40;
const PAIRS_PER_BATCH: usize = TILE_SIDE * TILE_SIDE;

/// A bounded unit of alignment work. Indices refer to the caller's protein or
/// candidate-pair collection; batches never copy sequences or expand all pairs.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
#[serde(
    tag = "kind",
    rename_all = "camelCase",
    rename_all_fields = "camelCase"
)]
pub enum ProteinBatch {
    Cartesian {
        query: Range<usize>,
        target: Range<usize>,
    },
    Selected {
        pairs: Range<usize>,
    },
}

impl ProteinBatch {
    pub fn pair_count(&self) -> usize {
        match self {
            Self::Cartesian { query, target } => query.len() * target.len(),
            Self::Selected { pairs } => pairs.len(),
        }
    }

    pub fn is_valid(&self, protein_count: usize, candidate_count: usize) -> bool {
        let valid = |range: &Range<usize>, end| range.start <= range.end && range.end <= end;
        match self {
            Self::Cartesian { query, target } => {
                valid(query, protein_count) && valid(target, protein_count)
            }
            Self::Selected { pairs } => valid(pairs, candidate_count),
        }
    }

    pub fn pairs<'a>(
        &self,
        candidates: &'a [ProteinPair],
    ) -> impl Iterator<Item = ProteinPair> + 'a {
        let (query, target, selected) = match self {
            Self::Cartesian { query, target } => (query.clone(), target.clone(), &[][..]),
            Self::Selected { pairs } => (0..0, 0..0, &candidates[pairs.clone()]),
        };
        query
            .flat_map(move |query_index| {
                target.clone().map(move |target_index| ProteinPair {
                    query_index,
                    target_index,
                })
            })
            .chain(selected.iter().copied())
    }
}

/// Plan cross-cluster work, including only pairs involving a newly added cluster.
/// For filtered inputs, the caller supplies the count of already selected pairs.
fn protein_batches(
    cluster_sizes: &[usize],
    first_new_cluster: usize,
    candidate_count: Option<usize>,
) -> Vec<ProteinBatch> {
    if let Some(count) = candidate_count {
        return (0..count)
            .step_by(PAIRS_PER_BATCH)
            .map(|start| ProteinBatch::Selected {
                pairs: start..(start + PAIRS_PER_BATCH).min(count),
            })
            .collect();
    }
    let mut offsets = vec![0];
    for size in cluster_sizes {
        offsets.push(offsets.last().unwrap() + size);
    }
    let mut batches = Vec::new();
    for query in 0..cluster_sizes.len() {
        for target in (query + 1).max(first_new_cluster)..cluster_sizes.len() {
            for query_start in (offsets[query]..offsets[query + 1]).step_by(TILE_SIDE) {
                for target_start in (offsets[target]..offsets[target + 1]).step_by(TILE_SIDE) {
                    batches.push(ProteinBatch::Cartesian {
                        query: query_start..(query_start + TILE_SIDE).min(offsets[query + 1]),
                        target: target_start..(target_start + TILE_SIDE).min(offsets[target + 1]),
                    });
                }
            }
        }
    }
    batches
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct AlignmentCounts {
    pub possible_pairs: usize,
    pub planned_pairs: usize,
    pub batch_count: usize,
}

/// Shared preparation for native execution and browser worker dispatch.
/// Sequences borrow the parsed clusters; only filtered candidate pairs are stored.
pub struct PreparedAlignment<'a> {
    pub proteins: Vec<&'a [u8]>,
    pub references: Vec<GeneRef>,
    pub candidates: Option<Vec<ProteinPair>>,
    pub batches: Vec<ProteinBatch>,
    pub counts: AlignmentCounts,
}

impl<'a> PreparedAlignment<'a> {
    pub fn new(
        clusters: &'a [Cluster],
        first_new_cluster: usize,
        options: AnalysisOptions,
        progress: impl FnMut(usize, usize),
    ) -> Result<Self, AnalysisOptionsError> {
        options.validate()?;
        let mut proteins = Vec::new();
        let mut references = Vec::new();
        let mut sizes = Vec::with_capacity(clusters.len());
        for (cluster, data) in clusters.iter().enumerate() {
            let start = proteins.len();
            for (locus, data) in data.loci.iter().enumerate() {
                for (gene, data) in data.genes.iter().enumerate() {
                    proteins.push(data.translation.as_bytes());
                    references.push(GeneRef {
                        cluster,
                        locus,
                        gene,
                    });
                }
            }
            sizes.push(proteins.len() - start);
        }
        let candidates = options.prefilter.map(|prefilter| {
            let membership = references.iter().map(|r| r.cluster).collect::<Vec<_>>();
            let mut pairs = cluster_kmer_candidate_pairs_with_progress(
                &proteins,
                &membership,
                options.identity_cutoff,
                prefilter,
                progress,
            );
            pairs.retain(|pair| references[pair.target_index].cluster >= first_new_cluster);
            pairs
        });
        let batches = protein_batches(&sizes, first_new_cluster, candidates.as_ref().map(Vec::len));
        let mut preceding = 0;
        let mut possible_pairs = 0;
        for (cluster, size) in sizes.into_iter().enumerate() {
            if cluster >= first_new_cluster {
                possible_pairs += preceding * size;
            }
            preceding += size;
        }
        let counts = AlignmentCounts {
            possible_pairs,
            planned_pairs: candidates.as_ref().map_or(possible_pairs, Vec::len),
            batch_count: batches.len(),
        };
        Ok(Self {
            proteins,
            references,
            candidates,
            batches,
            counts,
        })
    }

    pub(crate) fn run(
        &self,
        batch: &ProteinBatch,
        aligner: &mut ProteinBatchAligner,
        cutoff: f32,
    ) -> Result<(Vec<Link>, usize), AnalysisOptionsError> {
        let result = aligner.analyse_pairs(
            &self.proteins,
            &self.proteins,
            batch.pairs(self.candidates.as_deref().unwrap_or_default()),
            cutoff,
            |_| {},
        )?;
        let links = result
            .matches
            .into_iter()
            .map(|matched| Link {
                query: self.references[matched.query_index],
                target: self.references[matched.target_index],
                identity: matched.identity,
                similarity: matched.similarity,
            })
            .collect();
        Ok((links, result.aligned_pairs))
    }
}

#[derive(Debug, Default)]
pub struct ProteinBatchResult {
    pub matches: Vec<ProteinTileMatch>,
    pub aligned_pairs: usize,
}

/// Reusable alignment buffers for a native thread or browser worker.
#[derive(Default)]
pub struct ProteinBatchAligner {
    aligner: Option<ProteinAligner>,
    query_capacity: usize,
    target_capacity: usize,
}

impl ProteinBatchAligner {
    /// Compare an explicit or lazily generated sequence of pair indices. Invalid
    /// indices and pairs below the exact length bound are skipped. Progress counts
    /// examined pairs, including skipped pairs, every 100 pairs and at completion.
    pub fn analyse_pairs(
        &mut self,
        query: &[&[u8]],
        target: &[&[u8]],
        pairs: impl IntoIterator<Item = ProteinPair>,
        identity_cutoff: f32,
        mut progress: impl FnMut(usize),
    ) -> Result<ProteinBatchResult, AnalysisOptionsError> {
        AnalysisOptions {
            identity_cutoff,
            prefilter: None,
        }
        .validate()?;
        let mut result = ProteinBatchResult::default();
        let mut processed = 0;
        for pair in pairs {
            processed += 1;
            if processed % 100 == 0 {
                progress(processed);
            }
            let (Some(query), Some(target)) =
                (query.get(pair.query_index), target.get(pair.target_index))
            else {
                continue;
            };
            if !can_reach_identity(query.len(), target.len(), identity_cutoff) {
                continue;
            }
            if self.aligner.is_none()
                || query.len() > self.query_capacity
                || target.len() > self.target_capacity
            {
                self.query_capacity = self.query_capacity.max(query.len());
                self.target_capacity = self.target_capacity.max(target.len());
                self.aligner = Some(ProteinAligner::with_capacity(
                    self.query_capacity,
                    self.target_capacity,
                ));
            }
            result.aligned_pairs += 1;
            let matched = self.aligner.as_mut().unwrap().compare(query, target);
            if matched.identity >= identity_cutoff {
                result.matches.push(ProteinTileMatch {
                    query_index: pair.query_index,
                    target_index: pair.target_index,
                    identity: matched.identity,
                    similarity: matched.similarity,
                });
            }
        }
        if processed % 100 != 0 {
            progress(processed);
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn batches_bound_work_without_omitting_or_repeating_pairs() {
        let sizes = [41, 0, 83, 2];
        let mut membership = Vec::new();
        for (cluster, &size) in sizes.iter().enumerate() {
            membership.extend(std::iter::repeat_n(cluster, size));
        }
        for first_new in [0, 3, 4] {
            let batches = protein_batches(&sizes, first_new, None);
            assert!(batches.iter().all(|b| b.pair_count() <= PAIRS_PER_BATCH));
            let mut actual = batches
                .iter()
                .flat_map(|batch| batch.pairs(&[]))
                .map(|p| (p.query_index, p.target_index))
                .collect::<Vec<_>>();
            actual.sort_unstable();
            let expected = (0..membership.len())
                .flat_map(|q| {
                    let membership = &membership;
                    (q + 1..membership.len())
                        .filter(move |&t| {
                            membership[q] < membership[t] && membership[t] >= first_new
                        })
                        .map(move |t| (q, t))
                })
                .collect::<Vec<_>>();
            assert_eq!(actual, expected);
        }
        let pairs = vec![
            ProteinPair {
                query_index: 0,
                target_index: 1
            };
            3201
        ];
        let batches = protein_batches(&sizes, 0, Some(pairs.len()));
        assert_eq!(
            batches
                .iter()
                .map(ProteinBatch::pair_count)
                .collect::<Vec<_>>(),
            [1600, 1600, 1]
        );
        assert_eq!(
            batches
                .iter()
                .flat_map(|b| b.pairs(&pairs))
                .collect::<Vec<_>>(),
            pairs
        );
    }

    #[test]
    fn reusable_buffers_grow_across_batches_and_count_skipped_pairs() {
        let mut aligner = ProteinBatchAligner::default();
        assert!(
            aligner
                .analyse_pairs(&[b"MST"], &[b"MST"], [], f32::NAN, |_| panic!(
                    "invalid cutoff must not start alignment"
                ))
                .is_err()
        );
        for protein in [b"MST".as_slice(), b"MSTAVKLMDEFGH", b"MS"] {
            let result = aligner
                .analyse_pairs(
                    &[protein],
                    &[protein],
                    [ProteinPair {
                        query_index: 0,
                        target_index: 0,
                    }],
                    0.3,
                    |_| {},
                )
                .unwrap();
            assert_eq!(result.aligned_pairs, 1);
            assert_eq!(result.matches[0].identity, 1.0);
        }
        let mut progress = Vec::new();
        let result = aligner
            .analyse_pairs(
                &[b"MST"],
                &[b"MSTAVKLMDEFGH"],
                [
                    ProteinPair {
                        query_index: 0,
                        target_index: 0,
                    },
                    ProteinPair {
                        query_index: 9,
                        target_index: 0,
                    },
                ],
                0.9,
                |count| progress.push(count),
            )
            .unwrap();
        assert_eq!(result.aligned_pairs, 0);
        assert!(result.matches.is_empty());
        assert_eq!(progress, [2]);
    }
}
