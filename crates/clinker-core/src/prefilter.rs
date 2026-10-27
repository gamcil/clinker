use std::{
    collections::{BTreeMap, HashMap},
    hash::Hash,
};

use crate::{ProteinPair, can_reach_identity, protein::normalise_protein};

/// Approximate selection by distinct protein words; candidates still need alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq, serde::Deserialize, serde::Serialize)]
#[serde(rename_all = "camelCase", deny_unknown_fields)]
pub struct KmerPrefilter {
    pub kmer_size: usize,
    pub min_shared_kmers: usize,
}

impl Default for KmerPrefilter {
    fn default() -> Self {
        Self {
            kmer_size: 3,
            min_shared_kmers: 3,
        }
    }
}

#[cfg(test)]
pub(crate) fn kmer_candidate_pairs(
    query: &[&[u8]],
    target: &[&[u8]],
    identity_cutoff: f32,
    prefilter: KmerPrefilter,
) -> Vec<ProteinPair> {
    kmer_candidate_pairs_with_progress(query, target, identity_cutoff, prefilter, |_| {})
}

/// Report searched query-target pairs, including those rejected by the filter.
#[cfg(test)]
pub(crate) fn kmer_candidate_pairs_with_progress(
    query: &[&[u8]],
    target: &[&[u8]],
    identity_cutoff: f32,
    prefilter: KmerPrefilter,
    mut progress: impl FnMut(usize),
) -> Vec<ProteinPair> {
    let proteins = query.iter().chain(target).copied().collect::<Vec<_>>();
    let clusters = (0..proteins.len())
        .map(|i| usize::from(i >= query.len()))
        .collect::<Vec<_>>();
    cluster_kmer_candidate_pairs_with_progress(
        &proteins,
        &clusters,
        identity_cutoff,
        prefilter,
        |completed, total| {
            if completed > 0 || total == 0 {
                progress(completed);
            }
        },
    )
    .into_iter()
    .map(|pair| ProteinPair {
        query_index: pair.query_index,
        target_index: pair.target_index - query.len(),
    })
    .collect()
}

/// Select cross-cluster pairs, caching each protein's words and indexing one target
/// cluster at a time. Output indices refer to `proteins`, sorted by cluster pair
/// then protein indices, with the lower cluster ID as query. `clusters` must
/// contain one ID per protein.
pub(crate) fn cluster_kmer_candidate_pairs_with_progress(
    proteins: &[&[u8]],
    clusters: &[usize],
    identity_cutoff: f32,
    prefilter: KmerPrefilter,
    progress: impl FnMut(usize, usize),
) -> Vec<ProteinPair> {
    // Normalised residues are letters A–Z: five bits each, with no hashing.
    // Six residues fit in u32; seven and eight require u64.
    if prefilter.kmer_size <= 6 {
        cluster_candidates::<u32>(proteins, clusters, identity_cutoff, prefilter, progress)
    } else {
        cluster_candidates::<u64>(proteins, clusters, identity_cutoff, prefilter, progress)
    }
}

fn cluster_candidates<K: Copy + Ord + Hash + TryFrom<u64>>(
    proteins: &[&[u8]],
    clusters: &[usize],
    identity_cutoff: f32,
    prefilter: KmerPrefilter,
    mut progress: impl FnMut(usize, usize),
) -> Vec<ProteinPair> {
    assert_eq!(proteins.len(), clusters.len());
    let mut by_cluster = BTreeMap::<usize, Vec<usize>>::new();
    for (index, &cluster) in clusters.iter().enumerate() {
        by_cluster.entry(cluster).or_default().push(index);
    }
    let groups = by_cluster.into_values().collect::<Vec<_>>();
    let mut preceding = 0;
    let total = groups
        .iter()
        .map(|group| {
            let pairs = preceding * group.len();
            preceding += group.len();
            pairs
        })
        .sum();
    progress(0, total);
    let words = proteins
        .iter()
        .map(|protein| distinct_kmers::<K>(protein, prefilter.kmer_size))
        .collect::<Vec<_>>();
    let mut candidates = Vec::new();
    let mut counts = vec![0; proteins.len()];
    let mut touched = Vec::new();
    let mut completed = 0;
    for (target_group, targets) in groups.iter().enumerate().skip(1) {
        let mut index = HashMap::<K, Vec<usize>>::new();
        for &target in targets {
            for &word in &words[target] {
                index.entry(word).or_default().push(target);
            }
        }
        for queries in &groups[..target_group] {
            for &query in queries {
                count_shared_words(&words[query], &index, &mut counts, &mut touched);
                for target in touched.drain(..) {
                    let shared = std::mem::take(&mut counts[target]);
                    if shared >= prefilter.min_shared_kmers
                        && can_reach_identity(
                            proteins[query].len(),
                            proteins[target].len(),
                            identity_cutoff,
                        )
                    {
                        candidates.push(ProteinPair {
                            query_index: query,
                            target_index: target,
                        });
                    }
                }
                completed += targets.len();
                progress(completed, total);
            }
        }
    }
    candidates.sort_unstable_by_key(|pair| {
        (
            clusters[pair.query_index],
            clusters[pair.target_index],
            pair.query_index,
            pair.target_index,
        )
    });
    candidates
}

// A separate hot loop lets WASM engines optimise repeated queries during a run.
#[cfg_attr(target_arch = "wasm32", inline(never))]
fn count_shared_words<K: Eq + Hash>(
    words: &[K],
    index: &HashMap<K, Vec<usize>>,
    counts: &mut [usize],
    touched: &mut Vec<usize>,
) {
    for word in words {
        if let Some(targets) = index.get(word) {
            for &target in targets {
                if counts[target] == 0 {
                    touched.push(target);
                }
                counts[target] += 1;
            }
        }
    }
}

fn distinct_kmers<K: Ord + TryFrom<u64>>(protein: &[u8], kmer_size: usize) -> Box<[K]> {
    let mut normalised = Vec::with_capacity(protein.len());
    normalise_protein(protein, &mut normalised);
    let mut words = normalised
        .windows(kmer_size)
        .map(|word| {
            let key = word.iter().fold(0_u64, |key, &amino_acid| {
                (key << 5) | u64::from(amino_acid - b'A')
            });
            K::try_from(key)
                .unwrap_or_else(|_| unreachable!("word size selects a sufficient key width"))
        })
        .collect::<Vec<_>>();
    words.sort_unstable();
    words.dedup();
    words.into_boxed_slice()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cached_selection_matches_independent_word_intersections() {
        let proteins = [
            b"MSTAVKDLGYRQ".as_slice(),
            b"mstAVKDLGYRQ",
            b"OOOOOUUUUUXX",
            b"XXXXXXXXXXXX",
            b"MSTAVKD",
            b"ECCCCCCC",
            b"ACCCCCCC",
        ];
        let clusters = [4, 1, 4, 7, 1, 7, 1];
        for kmer_size in [1, 3, 5, 6, 7, 8] {
            let settings = KmerPrefilter {
                kmer_size,
                min_shared_kmers: 2,
            };
            let mut expected = Vec::new();
            for query in 0..proteins.len() {
                for target in 0..proteins.len() {
                    if clusters[query] >= clusters[target]
                        || !can_reach_identity(proteins[query].len(), proteins[target].len(), 0.8)
                    {
                        continue;
                    }
                    let byte_words = |protein: &[u8]| {
                        let mut normalised = Vec::new();
                        normalise_protein(protein, &mut normalised);
                        normalised
                            .windows(kmer_size)
                            .map(<[u8]>::to_vec)
                            .collect::<std::collections::BTreeSet<_>>()
                    };
                    let shared = byte_words(proteins[query])
                        .intersection(&byte_words(proteins[target]))
                        .count();
                    if shared >= settings.min_shared_kmers {
                        expected.push(ProteinPair {
                            query_index: query,
                            target_index: target,
                        });
                    }
                }
            }
            expected.sort_unstable_by_key(|pair| {
                (
                    clusters[pair.query_index],
                    clusters[pair.target_index],
                    pair.query_index,
                    pair.target_index,
                )
            });
            let mut final_progress = (0, 0);
            let actual = cluster_kmer_candidate_pairs_with_progress(
                &proteins,
                &clusters,
                0.8,
                settings,
                |done, total| final_progress = (done, total),
            );
            assert_eq!(actual, expected);
            assert_eq!(final_progress, (16, 16));
        }
    }
}
