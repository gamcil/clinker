//! Automatic locus ordering and orientation for plot output.

use std::collections::{BTreeMap, BTreeSet};

use crate::{Analysis, GeneRef, Link};

#[derive(Debug, Clone, Copy)]
pub(crate) struct LocusPlacement {
    pub locus: usize,
    pub reversed: bool,
    /// Absolute display offset in locus coordinates. `None` is used only while
    /// an unanchored locus is waiting for deterministic fallback placement.
    pub offset: Option<f64>,
    anchor: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
struct PositionedGene {
    reference: GeneRef,
    position: f64,
}

/// Arrange each cluster's loci against the already placed clusters.
#[cfg(test)]
pub(crate) fn auto_locus_layouts(
    analysis: &Analysis,
    cluster_order: &[usize],
) -> Vec<Vec<LocusPlacement>> {
    auto_locus_layouts_with_progress(analysis, cluster_order, |_| {})
}

pub(crate) fn auto_locus_layouts_with_progress(
    analysis: &Analysis,
    cluster_order: &[usize],
    mut progress: impl FnMut(usize),
) -> Vec<Vec<LocusPlacement>> {
    let scores = link_scores(&analysis.links);
    let neighbours = link_neighbours(&scores);
    let mut placed = Vec::<Vec<PositionedGene>>::new();
    let mut placed_positions = BTreeMap::<GeneRef, f64>::new();
    let mut placed_loci = BTreeMap::<(usize, usize), usize>::new();
    let mut alignment_row = Vec::new();
    let mut layouts = Vec::with_capacity(cluster_order.len());
    for &cluster_index in cluster_order {
        let cluster = &analysis.clusters[cluster_index];
        let mut layout = cluster
            .loci
            .iter()
            .enumerate()
            .map(|(locus, locus_data)| {
                let forward = positioned_genes(cluster_index, locus, locus_data, false);
                let reverse = positioned_genes(cluster_index, locus, locus_data, true);
                let references = forward
                    .iter()
                    .filter_map(|gene| neighbours.get(&gene.reference))
                    .flatten()
                    .filter_map(|(gene, _)| placed_loci.get(&(gene.cluster, gene.locus)).copied())
                    .collect::<BTreeSet<_>>();
                let mut best_score = |target: &[PositionedGene]| {
                    references
                        .iter()
                        .map(|&index| {
                            alignment_score(target, &placed[index], &scores, &mut alignment_row)
                        })
                        .fold(0.0_f64, f64::max)
                };
                let forward_score = best_score(&forward);
                let reverse_score = best_score(&reverse);
                let reversed = reverse_score > forward_score;
                let positioned = if reversed { &reverse } else { &forward };
                let placement = consensus_placement(positioned, &neighbours, &placed_positions);
                LocusPlacement {
                    locus,
                    reversed,
                    offset: placement.map(|(_, offset)| offset),
                    anchor: placement.map(|(anchor, _)| anchor),
                }
            })
            .collect::<Vec<_>>();
        layout.sort_by(|left, right| match (left.anchor, right.anchor) {
            (Some(left_anchor), Some(right_anchor)) => left_anchor
                .total_cmp(&right_anchor)
                .then(left.locus.cmp(&right.locus)),
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => left.locus.cmp(&right.locus),
        });
        place_unanchored_loci(cluster, &mut layout);
        for placement in &layout {
            let locus = &cluster.loci[placement.locus];
            let offset = placement.offset.unwrap_or(0.0);
            let positioned =
                positioned_genes(cluster_index, placement.locus, locus, placement.reversed)
                    .into_iter()
                    .map(|mut gene| {
                        gene.position += offset;
                        gene
                    })
                    .collect::<Vec<_>>();
            placed_positions.extend(
                positioned
                    .iter()
                    .map(|gene| (gene.reference, gene.position)),
            );
            placed_loci.insert((cluster_index, placement.locus), placed.len());
            placed.push(positioned);
        }
        layouts.push(layout);
        progress(layouts.len());
    }
    center_disconnected_components(analysis, cluster_order, &mut layouts);
    if cluster_order.is_empty() {
        progress(0);
    }
    layouts
}

/// Link-disconnected components have no evidence-defined relative x-position.
/// Preserve the first (normally most connected) component and centre every
/// other component on it so arbitrary locus lengths cannot create large
/// horizontal jumps between otherwise coherent blocks.
fn center_disconnected_components(
    analysis: &Analysis,
    cluster_order: &[usize],
    layouts: &mut [Vec<LocusPlacement>],
) {
    let Some(&reference_cluster) = cluster_order.first() else {
        return;
    };
    let mut parents = (0..analysis.clusters.len()).collect::<Vec<_>>();
    for link in &analysis.links {
        union_components(&mut parents, link.query.cluster, link.target.cluster);
    }

    let mut extents = BTreeMap::<usize, (f64, f64)>::new();
    for (&cluster_index, layout) in cluster_order.iter().zip(layouts.iter()) {
        let root = component_root(&parents, cluster_index);
        for placement in layout {
            let Some(offset) = placement.offset else {
                continue;
            };
            let locus = &analysis.clusters[cluster_index].loci[placement.locus];
            let end = offset + (locus.end - locus.start) as f64;
            extents
                .entry(root)
                .and_modify(|extent| {
                    extent.0 = extent.0.min(offset);
                    extent.1 = extent.1.max(end);
                })
                .or_insert((offset, end));
        }
    }

    let reference_root = component_root(&parents, reference_cluster);
    let Some(&(reference_start, reference_end)) = extents.get(&reference_root) else {
        return;
    };
    let reference_centre = (reference_start + reference_end) / 2.0;
    let shifts = extents
        .iter()
        .map(|(&root, &(start, end))| (root, reference_centre - (start + end) / 2.0))
        .collect::<BTreeMap<_, _>>();

    for (&cluster_index, layout) in cluster_order.iter().zip(layouts.iter_mut()) {
        let shift = shifts[&component_root(&parents, cluster_index)];
        for placement in layout {
            if let Some(offset) = &mut placement.offset {
                *offset += shift;
            }
        }
    }
}

fn component_root(parents: &[usize], mut cluster: usize) -> usize {
    while parents[cluster] != cluster {
        cluster = parents[cluster];
    }
    cluster
}

fn union_components(parents: &mut [usize], left: usize, right: usize) {
    let left_root = component_root(parents, left);
    let right_root = component_root(parents, right);
    if left_root != right_root {
        parents[right_root] = left_root;
    }
}

/// Give loci without a homology anchor a real position before they can become
/// references for later clusters. Leaving these positions to the renderer
/// would make the displayed packing differ from the zero offset assumed by
/// the following dynamic-programming pass.
fn place_unanchored_loci(cluster: &crate::Cluster, layout: &mut [LocusPlacement]) {
    let gap = typical_gene_span(cluster);
    let rightmost_anchor = layout
        .iter()
        .filter_map(|placement| {
            placement.offset.map(|offset| {
                offset
                    + (cluster.loci[placement.locus].end - cluster.loci[placement.locus].start)
                        as f64
            })
        })
        .max_by(f64::total_cmp);
    let mut cursor = rightmost_anchor.map_or(0.0, |right| right + gap);

    for placement in layout
        .iter_mut()
        .filter(|placement| placement.offset.is_none())
    {
        placement.offset = Some(cursor);
        cursor +=
            (cluster.loci[placement.locus].end - cluster.loci[placement.locus].start) as f64 + gap;
    }
}

/// Use roughly one gene of separation so fallback packing scales naturally
/// across compact and large-coordinate records without depending on a
/// renderer-specific pixel spacing.
fn typical_gene_span(cluster: &crate::Cluster) -> f64 {
    let mut spans = cluster
        .loci
        .iter()
        .flat_map(|locus| &locus.genes)
        .filter_map(|gene| gene.end.checked_sub(gene.start).filter(|span| *span > 0))
        .collect::<Vec<_>>();
    spans.sort_unstable();
    spans.get(spans.len() / 2).copied().unwrap_or(1) as f64
}

fn positioned_genes(
    cluster: usize,
    locus: usize,
    locus_data: &crate::Locus,
    reversed: bool,
) -> Vec<PositionedGene> {
    // The renderer normalizes every locus to a zero origin, so offsets must
    // use the same locus-local coordinate system rather than GenBank's raw
    // record coordinates.
    let centre = |gene: &crate::Gene| {
        (gene.start - locus_data.start) as f64 + (gene.end - gene.start) as f64 / 2.0
    };
    let display_position = |gene: &crate::Gene| {
        if reversed {
            (locus_data.end - locus_data.start) as f64 - centre(gene)
        } else {
            centre(gene)
        }
    };
    (0..locus_data.genes.len())
        .map(|index| {
            if reversed {
                locus_data.genes.len() - 1 - index
            } else {
                index
            }
        })
        .map(|gene| PositionedGene {
            reference: GeneRef {
                cluster,
                locus,
                gene,
            },
            position: display_position(&locus_data.genes[gene]),
        })
        .collect()
}

fn link_scores(links: &[Link]) -> BTreeMap<(GeneRef, GeneRef), f64> {
    let mut scores = BTreeMap::new();
    for link in links {
        let key = ordered_pair(link.query, link.target);
        scores
            .entry(key)
            .and_modify(|score: &mut f64| *score = (*score).max(link.identity as f64))
            .or_insert(link.identity as f64);
    }
    scores
}

fn link_neighbours(
    scores: &BTreeMap<(GeneRef, GeneRef), f64>,
) -> BTreeMap<GeneRef, Vec<(GeneRef, f64)>> {
    let mut neighbours = BTreeMap::<GeneRef, Vec<(GeneRef, f64)>>::new();
    for (&(left, right), &score) in scores {
        neighbours.entry(left).or_default().push((right, score));
        neighbours.entry(right).or_default().push((left, score));
    }
    neighbours
}

fn ordered_pair(one: GeneRef, two: GeneRef) -> (GeneRef, GeneRef) {
    if one <= two { (one, two) } else { (two, one) }
}

const MISMATCH_PENALTY: f64 = -1.0;
const GAP_PENALTY: f64 = -0.25;

fn consensus_placement(
    target: &[PositionedGene],
    neighbours: &BTreeMap<GeneRef, Vec<(GeneRef, f64)>>,
    placed_positions: &BTreeMap<GeneRef, f64>,
) -> Option<(f64, f64)> {
    let mut anchors = Vec::new();
    let mut offsets = Vec::new();
    for target_gene in target {
        let Some(linked_genes) = neighbours.get(&target_gene.reference) else {
            continue;
        };
        for &(reference, weight) in linked_genes {
            let Some(&reference_position) = placed_positions.get(&reference) else {
                continue;
            };
            anchors.push((reference_position, weight));
            offsets.push((reference_position - target_gene.position, weight));
        }
    }
    Some((weighted_median(anchors)?, weighted_median(offsets)?))
}

fn weighted_median(mut values: Vec<(f64, f64)>) -> Option<f64> {
    values.sort_by(|left, right| left.0.total_cmp(&right.0));
    let midpoint = values.iter().map(|(_, weight)| weight).sum::<f64>() / 2.0;
    let mut cumulative = 0.0;
    for (value, weight) in &values {
        cumulative += weight;
        if cumulative >= midpoint {
            return Some(*value);
        }
    }
    values.last().map(|(value, _)| *value)
}

// Orientation needs only the best collinear-link score. Placement uses the
// consensus anchors, so no traceback or full score matrix is needed here.
fn alignment_score(
    target: &[PositionedGene],
    reference: &[PositionedGene],
    scores: &BTreeMap<(GeneRef, GeneRef), f64>,
    row: &mut Vec<f64>,
) -> f64 {
    row.resize(reference.len() + 1, 0.0);
    row.fill(0.0);
    let mut best = 0.0_f64;
    for target_gene in target {
        let mut diagonal = 0.0;
        for (j, reference_gene) in reference.iter().enumerate() {
            let above = row[j + 1];
            let match_score = scores
                .get(&ordered_pair(
                    target_gene.reference,
                    reference_gene.reference,
                ))
                .copied()
                .unwrap_or(MISMATCH_PENALTY);
            let value = (diagonal + match_score)
                .max(above + GAP_PENALTY)
                .max(row[j] + GAP_PENALTY)
                .max(0.0);
            row[j + 1] = value;
            diagonal = above;
            best = best.max(value);
        }
    }
    best
}

#[cfg(test)]
mod tests {
    use super::{auto_locus_layouts, weighted_median};
    use crate::{Analysis, Cluster, Gene, GeneRef, Link, Locus};

    fn cluster(loci: usize, genes: usize) -> Cluster {
        Cluster {
            name: "test".into(),
            source_file: None,
            organism: None,
            loci: (0..loci)
                .map(|locus| {
                    let start = locus as u64 * 100;
                    Locus {
                        name: "locus".into(),
                        start,
                        end: start + genes as u64,
                        genes: (0..genes)
                            .map(|gene| Gene {
                                label: gene.to_string(),
                                names: Vec::new(),
                                start: start + gene as u64,
                                end: start + gene as u64 + 1,
                                strand: 1,
                                translation: "M".into(),
                            })
                            .collect(),
                    }
                })
                .collect(),
        }
    }

    fn link(query: GeneRef, target: GeneRef) -> Link {
        Link {
            query,
            target,
            identity: 1.0,
            similarity: 1.0,
        }
    }

    #[test]
    fn orders_loci_by_anchor_and_flips_reversed_gene_order() {
        let analysis = Analysis {
            clusters: vec![cluster(1, 4), cluster(2, 2)],
            links: vec![
                link(
                    GeneRef {
                        cluster: 0,
                        locus: 0,
                        gene: 0,
                    },
                    GeneRef {
                        cluster: 1,
                        locus: 1,
                        gene: 1,
                    },
                ),
                link(
                    GeneRef {
                        cluster: 0,
                        locus: 0,
                        gene: 1,
                    },
                    GeneRef {
                        cluster: 1,
                        locus: 1,
                        gene: 0,
                    },
                ),
                link(
                    GeneRef {
                        cluster: 0,
                        locus: 0,
                        gene: 2,
                    },
                    GeneRef {
                        cluster: 1,
                        locus: 0,
                        gene: 0,
                    },
                ),
                link(
                    GeneRef {
                        cluster: 0,
                        locus: 0,
                        gene: 3,
                    },
                    GeneRef {
                        cluster: 1,
                        locus: 0,
                        gene: 1,
                    },
                ),
            ],
        };
        let layouts = auto_locus_layouts(&analysis, &[0, 1]);
        assert_eq!(layouts[1][0].locus, 1);
        assert!(layouts[1][0].reversed);
        assert_eq!(layouts[1][0].offset, Some(0.0));
        assert_eq!(layouts[1][1].locus, 0);
        assert_eq!(layouts[1][1].offset, Some(2.0));
    }

    #[test]
    fn disconnected_root_uses_the_same_position_for_later_alignments() {
        let analysis = Analysis {
            clusters: vec![cluster(2, 2), cluster(1, 2)],
            links: vec![link(
                GeneRef {
                    cluster: 0,
                    locus: 1,
                    gene: 0,
                },
                GeneRef {
                    cluster: 1,
                    locus: 0,
                    gene: 0,
                },
            )],
        };

        let layouts = auto_locus_layouts(&analysis, &[0, 1]);
        let root_offset = layouts[0]
            .iter()
            .find(|placement| placement.locus == 1)
            .and_then(|placement| placement.offset)
            .expect("a disconnected root locus should receive an explicit offset");
        let aligned_offset = layouts[1][0]
            .offset
            .expect("a linked locus should receive an explicit offset");

        assert_eq!(aligned_offset, root_offset);
    }

    #[test]
    fn centres_disconnected_components_on_the_reference_component() {
        let analysis = Analysis {
            clusters: vec![cluster(1, 2), cluster(1, 6)],
            links: Vec::new(),
        };

        let layouts = auto_locus_layouts(&analysis, &[0, 1]);
        let first_centre = layouts[0][0].offset.unwrap() + 1.0;
        let second_centre = layouts[1][0].offset.unwrap() + 3.0;

        assert_eq!(second_centre, first_centre);
    }

    #[test]
    fn does_not_chain_an_alignment_across_separate_reference_loci() {
        let analysis = Analysis {
            clusters: vec![cluster(2, 2), cluster(1, 2)],
            links: vec![
                link(
                    GeneRef {
                        cluster: 0,
                        locus: 0,
                        gene: 0,
                    },
                    GeneRef {
                        cluster: 1,
                        locus: 0,
                        gene: 0,
                    },
                ),
                link(
                    GeneRef {
                        cluster: 0,
                        locus: 1,
                        gene: 1,
                    },
                    GeneRef {
                        cluster: 1,
                        locus: 0,
                        gene: 1,
                    },
                ),
            ],
        };

        let layouts = auto_locus_layouts(&analysis, &[0, 1]);
        let reference_offsets = layouts[0]
            .iter()
            .map(|placement| placement.offset.unwrap())
            .collect::<Vec<_>>();
        let target_offset = layouts[1][0].offset.unwrap();

        assert!(reference_offsets.contains(&target_offset));
    }

    #[test]
    fn consensus_offset_resists_a_distant_low_support_anchor() {
        assert_eq!(
            weighted_median(vec![(10.0, 0.9), (11.0, 0.8), (250.0, 0.3)]),
            Some(11.0),
        );
    }

    #[test]
    fn tied_paralogue_anchors_keep_forward_orientation_independent_of_link_order() {
        let gene = |cluster, gene| GeneRef {
            cluster,
            locus: 0,
            gene,
        };
        let mut analysis = Analysis {
            clusters: vec![cluster(1, 1), cluster(1, 2)],
            links: vec![link(gene(0, 0), gene(1, 0)), link(gene(0, 0), gene(1, 1))],
        };
        let layout = auto_locus_layouts(&analysis, &[0, 1]);
        assert!(!layout[1][0].reversed);
        assert_eq!(layout[1][0].offset, Some(-1.0));
        analysis.links.reverse();
        let reordered = auto_locus_layouts(&analysis, &[0, 1]);
        assert_eq!(reordered[1][0].reversed, layout[1][0].reversed);
        assert_eq!(reordered[1][0].offset, layout[1][0].offset);
    }
}
