//! Cluster-pair synteny scores and their distance-matrix representation.

use std::collections::BTreeMap;

use crate::hierarchy::ward_leaf_order;
use crate::{Analysis, Cluster, GeneRef, Link, build_gene_groups};

/// Weight applied to shared adjacent homology-group pairs.
pub const DEFAULT_CONTIGUITY_WEIGHT: f64 = 0.5;

impl Analysis {
    /// Calculate the clinker synteny score matrix.
    /// For a pair of clusters this is the sum of retained-link identities plus
    /// `contiguity_weight` times the number of shared adjacent group pairs.
    pub fn synteny_matrix(&self, contiguity_weight: f64) -> Vec<Vec<f64>> {
        let count = self.clusters.len();
        let mut matrix = vec![vec![0.0; count]; count];
        for query in 0..count {
            for target in query + 1..count {
                let links = self
                    .links
                    .iter()
                    .filter(|link| link.query.cluster == query && link.target.cluster == target)
                    .cloned()
                    .collect::<Vec<_>>();
                let score = synteny_score(
                    &self.clusters[query],
                    &self.clusters[target],
                    &links,
                    contiguity_weight,
                );
                matrix[query][target] = score;
                matrix[target][query] = score;
            }
        }
        matrix
    }

    /// Create a normalized distance matrix suitable for hierarchical ordering.
    pub fn synteny_distance_matrix(&self, contiguity_weight: f64) -> Vec<Vec<f64>> {
        let mut matrix = self.synteny_matrix(contiguity_weight);
        let maximum = matrix
            .iter()
            .flat_map(|row| row.iter())
            .copied()
            .fold(0.0_f64, f64::max);
        if maximum > 0.0 {
            for row in &mut matrix {
                for value in row {
                    *value /= maximum;
                }
            }
        }
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

    /// Order clusters by Ward linkage over normalized synteny distance.
    pub fn cluster_order(&self, contiguity_weight: f64) -> Vec<usize> {
        let count = self.clusters.len();
        if count < 2 || self.links.is_empty() {
            return (0..count).collect();
        }
        let matrix = self.synteny_distance_matrix(contiguity_weight);
        let mut condensed = Vec::with_capacity(count * (count - 1) / 2);
        for row in 0..count {
            condensed.extend_from_slice(&matrix[row][row + 1..]);
        }
        // Python reverses SciPy's leaf sequence before plotting.
        let mut order = ward_leaf_order(&condensed, count);
        order.reverse();
        order
    }
}

fn synteny_score(query: &Cluster, target: &Cluster, links: &[Link], contiguity_weight: f64) -> f64 {
    let homology = links
        .iter()
        .map(|link| f64::from(link.identity))
        .sum::<f64>();
    if links.is_empty() {
        return 0.0;
    }
    let mut groups = BTreeMap::new();
    for (index, group) in build_gene_groups(links).into_iter().enumerate() {
        for gene in group.genes {
            groups.insert(gene, index);
        }
    }
    let query_pairs = adjacent_group_counts(query, 0, &groups);
    let target_pairs = adjacent_group_counts(target, 1, &groups);
    let contiguity = query_pairs
        .iter()
        .filter_map(|(pair, count)| target_pairs.get(pair).map(|other| count.min(other)))
        .sum::<usize>();
    homology + contiguity_weight * contiguity as f64
}

fn adjacent_group_counts(
    cluster: &Cluster,
    cluster_index: usize,
    groups: &BTreeMap<GeneRef, usize>,
) -> BTreeMap<(usize, usize), usize> {
    let mut counts = BTreeMap::new();
    for (locus_index, locus) in cluster.loci.iter().enumerate() {
        for gene_index in 0..locus.genes.len().saturating_sub(1) {
            let one = GeneRef {
                cluster: cluster_index,
                locus: locus_index,
                gene: gene_index,
            };
            let two = GeneRef {
                cluster: cluster_index,
                locus: locus_index,
                gene: gene_index + 1,
            };
            if let (Some(one_group), Some(two_group)) = (groups.get(&one), groups.get(&two)) {
                *counts.entry((*one_group, *two_group)).or_default() += 1;
            }
        }
    }
    counts
}

#[cfg(test)]
mod tests {
    use crate::{Analysis, Cluster, Gene, GeneRef, Link, Locus};

    fn cluster(name: &str, genes: usize) -> Cluster {
        Cluster {
            name: name.into(),
            loci: vec![Locus {
                name: name.into(),
                start: 0,
                end: genes,
                genes: (0..genes)
                    .map(|index| Gene {
                        label: format!("{name}-{index}"),
                        names: Vec::new(),
                        start: index,
                        end: index + 1,
                        strand: 1,
                        translation: "M".into(),
                    })
                    .collect(),
            }],
        }
    }

    fn link(query_gene: usize, target_gene: usize, identity: f32) -> Link {
        Link {
            query: GeneRef {
                cluster: 0,
                locus: 0,
                gene: query_gene,
            },
            target: GeneRef {
                cluster: 1,
                locus: 0,
                gene: target_gene,
            },
            identity,
            similarity: identity,
        }
    }

    #[test]
    fn includes_homology_and_shared_adjacency_in_synteny_scores() {
        let analysis = Analysis {
            clusters: vec![cluster("one", 2), cluster("two", 2)],
            links: vec![link(0, 0, 0.8), link(1, 1, 0.9)],
        };
        assert!((analysis.synteny_matrix(0.5)[0][1] - 2.2).abs() < 1e-6);
        assert_eq!(
            analysis.synteny_distance_matrix(0.5),
            vec![vec![0.0, 0.0], vec![0.0, 0.0]]
        );
    }
}
