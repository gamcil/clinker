//! Cluster-pair similarity scores and their distance-matrix representation.

use std::collections::BTreeMap;

use serde::Serialize;

use crate::hierarchy::ward_leaf_order;
use crate::{Analysis, Cluster, GeneRef, Link};

/// Directional best-hit coverages and their containment-aware similarity.
#[derive(Debug, Clone, Copy, Default, PartialEq, Serialize)]
#[serde(rename_all = "camelCase")]
pub struct ClusterPairSimilarity {
    pub query_coverage: f64,
    pub target_coverage: f64,
    pub similarity: f64,
}

impl Analysis {
    /// Calculate containment-aware, per-gene best-hit similarity.
    ///
    /// Each gene contributes at most its strongest retained identity, weighted
    /// by its translated protein length. The larger directional coverage is
    /// used so a well-matched smaller cluster remains similar to a larger
    /// cluster that contains it.
    pub fn synteny_matrix(&self) -> Vec<Vec<f64>> {
        self.cluster_similarity_matrix()
            .into_iter()
            .map(|row| row.into_iter().map(|pair| pair.similarity).collect())
            .collect()
    }

    /// Calculate directional best-hit coverages for every cluster pair.
    ///
    /// At `[query][target]`, `query_coverage` is query → target and
    /// `target_coverage` is target → query. The mirrored entry swaps them.
    pub fn cluster_similarity_matrix(&self) -> Vec<Vec<ClusterPairSimilarity>> {
        self.cluster_similarity_matrix_with_progress(|_| {})
    }

    /// As [`Analysis::cluster_similarity_matrix`], reporting completed
    /// unordered cluster pairs.
    pub fn cluster_similarity_matrix_with_progress(
        &self,
        mut progress: impl FnMut(usize),
    ) -> Vec<Vec<ClusterPairSimilarity>> {
        let count = self.clusters.len();
        let mut matrix = vec![vec![ClusterPairSimilarity::default(); count]; count];
        let mut links_by_pair = BTreeMap::<(usize, usize), Vec<&Link>>::new();
        for link in &self.links {
            links_by_pair
                .entry((link.query.cluster, link.target.cluster))
                .or_default()
                .push(link);
        }
        let mut processed = 0;
        for query in 0..count {
            for target in query + 1..count {
                let links = links_by_pair
                    .get(&(query, target))
                    .map(Vec::as_slice)
                    .unwrap_or_default();
                let pair =
                    containment_similarity(&self.clusters[query], &self.clusters[target], links);
                matrix[query][target] = pair;
                matrix[target][query] = ClusterPairSimilarity {
                    query_coverage: pair.target_coverage,
                    target_coverage: pair.query_coverage,
                    similarity: pair.similarity,
                };
                processed += 1;
                progress(processed);
            }
        }
        if processed == 0 {
            progress(0);
        }
        matrix
    }

    /// Create a distance matrix suitable for hierarchical ordering.
    pub fn synteny_distance_matrix(&self) -> Vec<Vec<f64>> {
        let mut matrix = self.synteny_matrix();
        for (row_index, row) in matrix.iter_mut().enumerate() {
            for (column_index, value) in row.iter_mut().enumerate() {
                *value = if row_index == column_index {
                    0.0
                } else {
                    1.0 - *value
                };
            }
        }
        matrix
    }

    /// Order clusters by Ward linkage over best-hit containment distance.
    ///
    /// A complete reversal has the same neighbour-distance score, so orient
    /// the final sequence with its strongest linked component first. Component
    /// strength is the total pairwise similarity within it, then its number of
    /// clusters as a tie-breaker.
    pub fn cluster_order(&self) -> Vec<usize> {
        let similarities = self.cluster_similarity_matrix();
        self.cluster_order_from_similarity_matrix(&similarities)
    }

    /// Order clusters using an already calculated similarity matrix.
    /// Reusing it avoids recalculating every cluster pair during browser
    /// post-processing, where the same matrix is also exported as CSV.
    pub fn cluster_order_from_similarity_matrix(
        &self,
        matrix: &[Vec<ClusterPairSimilarity>],
    ) -> Vec<usize> {
        let count = self.clusters.len();
        if count < 2 || self.links.is_empty() {
            return (0..count).collect();
        }
        assert_eq!(matrix.len(), count);
        assert!(matrix.iter().all(|row| row.len() == count));
        let similarities = matrix
            .iter()
            .map(|row| row.iter().map(|pair| pair.similarity).collect::<Vec<_>>())
            .collect::<Vec<_>>();
        let mut condensed = Vec::with_capacity(count * (count - 1) / 2);
        for row in 0..count {
            condensed.extend(
                similarities[row][row + 1..]
                    .iter()
                    .map(|similarity| 1.0 - similarity),
            );
        }
        let mut order = ward_leaf_order(&condensed, count);
        orient_strongest_component_first(&mut order, &similarities);
        order
    }
}

fn orient_strongest_component_first(order: &mut [usize], similarities: &[Vec<f64>]) {
    let components = linked_components(similarities);
    let strengths = component_strengths(similarities, &components);
    let first = components[order[0]];
    let last = components[*order.last().expect("an order contains a cluster")];
    if compare_component_strengths(strengths[last], strengths[first]).is_gt() {
        order.reverse();
    }
}

fn linked_components(similarities: &[Vec<f64>]) -> Vec<usize> {
    let count = similarities.len();
    let mut parent = (0..count).collect::<Vec<_>>();
    for query in 0..count {
        for target in query + 1..count {
            if similarities[query][target] > 0.0 {
                union_components(&mut parent, query, target);
            }
        }
    }
    for cluster in 0..count {
        parent[cluster] = find_component(&mut parent, cluster);
    }
    parent
}

fn find_component(parent: &mut [usize], cluster: usize) -> usize {
    if parent[cluster] != cluster {
        parent[cluster] = find_component(parent, parent[cluster]);
    }
    parent[cluster]
}

fn union_components(parent: &mut [usize], one: usize, two: usize) {
    let one_root = find_component(parent, one);
    let two_root = find_component(parent, two);
    if one_root != two_root {
        // Retaining the lower root makes the component IDs reproducible.
        parent[one_root.max(two_root)] = one_root.min(two_root);
    }
}

fn component_strengths(similarities: &[Vec<f64>], components: &[usize]) -> Vec<(f64, usize)> {
    let mut strengths = vec![(0.0, 0); similarities.len()];
    for (cluster, &component) in components.iter().enumerate() {
        strengths[component].1 += 1;
        for other in cluster + 1..similarities.len() {
            if components[other] == component {
                strengths[component].0 += similarities[cluster][other];
            }
        }
    }
    strengths
}

fn compare_component_strengths(left: (f64, usize), right: (f64, usize)) -> std::cmp::Ordering {
    left.0
        .total_cmp(&right.0)
        .then_with(|| left.1.cmp(&right.1))
}

fn containment_similarity(
    query: &Cluster,
    target: &Cluster,
    links: &[&Link],
) -> ClusterPairSimilarity {
    let mut query_hits = BTreeMap::<GeneRef, f64>::new();
    let mut target_hits = BTreeMap::<GeneRef, f64>::new();
    for link in links {
        query_hits
            .entry(link.query)
            .and_modify(|identity| *identity = identity.max(f64::from(link.identity)))
            .or_insert(f64::from(link.identity));
        target_hits
            .entry(link.target)
            .and_modify(|identity| *identity = identity.max(f64::from(link.identity)))
            .or_insert(f64::from(link.identity));
    }
    let protein_length = |cluster: &Cluster, gene: GeneRef| {
        cluster.loci[gene.locus].genes[gene.gene].translation.len()
    };
    let coverage = |hits: &BTreeMap<GeneRef, f64>, cluster: &Cluster| {
        let total_length = cluster
            .loci
            .iter()
            .flat_map(|locus| locus.genes.iter())
            .map(|gene| gene.translation.len())
            .sum::<usize>();
        if total_length == 0 {
            0.0
        } else {
            hits.iter()
                .map(|(gene, identity)| *identity * protein_length(cluster, *gene) as f64)
                .sum::<f64>()
                / total_length as f64
        }
    };
    let query_coverage = coverage(&query_hits, query);
    let target_coverage = coverage(&target_hits, target);
    ClusterPairSimilarity {
        query_coverage,
        target_coverage,
        similarity: query_coverage.max(target_coverage),
    }
}

#[cfg(test)]
mod tests {
    use super::orient_strongest_component_first;
    use crate::{Analysis, Cluster, Gene, GeneRef, Link, Locus};

    fn cluster(name: &str, genes: usize) -> Cluster {
        cluster_with_lengths(name, &vec![1; genes])
    }

    fn cluster_with_lengths(name: &str, lengths: &[usize]) -> Cluster {
        Cluster {
            name: name.into(),
            loci: vec![Locus {
                name: name.into(),
                start: 0,
                end: lengths.iter().sum(),
                genes: lengths
                    .iter()
                    .enumerate()
                    .map(|(index, length)| Gene {
                        label: format!("{name}-{index}"),
                        names: Vec::new(),
                        start: index,
                        end: index + 1,
                        strand: 1,
                        translation: "M".repeat(*length),
                    })
                    .collect(),
            }],
        }
    }

    fn link(query_gene: usize, target_gene: usize, identity: f32) -> Link {
        link_between(0, query_gene, 1, target_gene, identity)
    }

    fn link_between(
        query_cluster: usize,
        query_gene: usize,
        target_cluster: usize,
        target_gene: usize,
        identity: f32,
    ) -> Link {
        Link {
            query: GeneRef {
                cluster: query_cluster,
                locus: 0,
                gene: query_gene,
            },
            target: GeneRef {
                cluster: target_cluster,
                locus: 0,
                gene: target_gene,
            },
            identity,
            similarity: identity,
        }
    }

    #[test]
    fn uses_the_best_hit_once_per_gene_and_normalizes_by_cluster_size() {
        let analysis = Analysis {
            clusters: vec![cluster("one", 2), cluster("two", 2)],
            links: vec![link(0, 0, 0.8), link(0, 1, 0.9), link(1, 1, 0.7)],
        };
        assert!((analysis.synteny_matrix()[0][1] - 0.85).abs() < 1e-6);
        assert!((analysis.synteny_distance_matrix()[0][1] - 0.15).abs() < 1e-6);
    }

    #[test]
    fn scores_cluster_pairs_after_the_first_input_pair() {
        let analysis = Analysis {
            clusters: vec![cluster("zero", 2), cluster("one", 2), cluster("two", 2)],
            links: vec![link_between(1, 0, 2, 0, 0.8), link_between(1, 1, 2, 1, 0.9)],
        };

        assert!((analysis.synteny_matrix()[1][2] - 0.85).abs() < 1e-6);
    }

    #[test]
    fn similarity_matrix_reports_each_cluster_pair() {
        let analysis = Analysis {
            clusters: vec![cluster("zero", 1), cluster("one", 1), cluster("two", 1)],
            links: Vec::new(),
        };
        let mut reports = Vec::new();

        let _ =
            analysis.cluster_similarity_matrix_with_progress(|processed| reports.push(processed));

        assert_eq!(reports, vec![1, 2, 3]);
    }

    #[test]
    fn weights_best_hits_by_protein_length() {
        let analysis = Analysis {
            clusters: vec![
                cluster_with_lengths("one", &[100, 10]),
                cluster_with_lengths("two", &[100, 10]),
            ],
            links: vec![link(0, 0, 0.9)],
        };

        assert!((analysis.synteny_matrix()[0][1] - 90.0 / 110.0).abs() < 1e-6);
    }

    #[test]
    fn puts_the_most_connected_component_at_the_top() {
        // Clusters 0–2 are the dense main block; 3–4 are a smaller pair.
        let similarities = vec![
            vec![0.0, 0.9, 0.8, 0.0, 0.0],
            vec![0.9, 0.0, 0.9, 0.0, 0.0],
            vec![0.8, 0.9, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0, 0.7],
            vec![0.0, 0.0, 0.0, 0.7, 0.0],
        ];
        let mut order = vec![3, 4, 2, 1, 0];

        orient_strongest_component_first(&mut order, &similarities);

        assert_eq!(order, vec![0, 1, 2, 4, 3]);
    }
}
