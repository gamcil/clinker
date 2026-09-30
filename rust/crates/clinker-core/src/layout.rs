//! Automatic locus ordering and orientation for plot output.

use std::collections::BTreeMap;

use crate::{Analysis, GeneRef, Link};

#[derive(Debug, Clone, Copy)]
pub(crate) struct LocusPlacement {
    pub locus: usize,
    pub reversed: bool,
    /// Absolute display offset for an anchored locus, in locus coordinates.
    /// `None` leaves an unlinked locus to the renderer's normal packing.
    pub offset: Option<f64>,
    anchor: Option<f64>,
}

#[derive(Debug, Clone, Copy)]
struct PositionedGene {
    reference: GeneRef,
    position: f64,
}

/// Arrange each cluster's loci against the already placed clusters.
pub(crate) fn auto_locus_layouts(
    analysis: &Analysis,
    cluster_order: &[usize],
) -> Vec<Vec<LocusPlacement>> {
    let scores = link_scores(&analysis.links);
    let mut placed = Vec::<PositionedGene>::new();
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
                let forward_alignment = align_locus(&forward, &placed, &scores);
                let reverse_alignment = align_locus(&reverse, &placed, &scores);
                let (reversed, alignment) = if reverse_alignment.score > forward_alignment.score {
                    (true, reverse_alignment)
                } else {
                    (false, forward_alignment)
                };
                LocusPlacement {
                    locus,
                    reversed,
                    offset: alignment.offset,
                    anchor: alignment.anchor,
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
        for placement in &layout {
            let locus = &cluster.loci[placement.locus];
            let offset = placement.offset.unwrap_or(0.0);
            placed.extend(
                positioned_genes(cluster_index, placement.locus, locus, placement.reversed)
                    .into_iter()
                    .map(|mut gene| {
                        gene.position += offset;
                        gene
                    }),
            );
        }
        layouts.push(layout);
    }
    layouts
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
    let centre =
        |gene: &crate::Gene| (gene.start + gene.end) as f64 / 2.0 - locus_data.start as f64;
    let display_position = |gene: &crate::Gene| {
        if reversed {
            (locus_data.end - locus_data.start) as f64 - centre(gene)
        } else {
            centre(gene)
        }
    };
    let indexes: Box<dyn Iterator<Item = usize>> = if reversed {
        Box::new((0..locus_data.genes.len()).rev())
    } else {
        Box::new(0..locus_data.genes.len())
    };
    indexes
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

fn ordered_pair(one: GeneRef, two: GeneRef) -> (GeneRef, GeneRef) {
    if one <= two { (one, two) } else { (two, one) }
}

struct Alignment {
    score: f64,
    anchor: Option<f64>,
    offset: Option<f64>,
}

const MISMATCH_PENALTY: f64 = -1.0;
const GAP_PENALTY: f64 = -0.25;

fn align_locus(
    target: &[PositionedGene],
    reference: &[PositionedGene],
    scores: &BTreeMap<(GeneRef, GeneRef), f64>,
) -> Alignment {
    if target.is_empty() || reference.is_empty() {
        return Alignment {
            score: 0.0,
            anchor: None,
            offset: None,
        };
    }
    let width = reference.len() + 1;
    let mut values = vec![0.0_f64; (target.len() + 1) * width];
    let mut trace = vec![0_u8; values.len()];
    for i in 1..=target.len() {
        for j in 1..=reference.len() {
            let match_score = scores
                .get(&ordered_pair(
                    target[i - 1].reference,
                    reference[j - 1].reference,
                ))
                .copied();
            let diagonal =
                values[(i - 1) * width + j - 1] + match_score.unwrap_or(MISMATCH_PENALTY);
            let up = values[(i - 1) * width + j] + GAP_PENALTY;
            let left = values[i * width + j - 1] + GAP_PENALTY;
            let (value, direction) = if diagonal > 0.0 && diagonal >= up && diagonal >= left {
                (diagonal, 1)
            } else if up > 0.0 && up >= left {
                (up, 2)
            } else if left > 0.0 {
                (left, 3)
            } else {
                (0.0, 0)
            };
            values[i * width + j] = value;
            trace[i * width + j] = direction;
        }
    }
    let (mut index, score) = values
        .iter()
        .enumerate()
        .max_by(|(_, a), (_, b)| a.total_cmp(b))
        .map(|(i, s)| (i, *s))
        .unwrap();
    let mut anchors = Vec::new();
    while trace[index] != 0 {
        let i = index / width;
        let j = index % width;
        match trace[index] {
            1 => {
                if let Some(weight) = scores.get(&ordered_pair(
                    target[i - 1].reference,
                    reference[j - 1].reference,
                )) {
                    anchors.push((
                        reference[j - 1].position,
                        reference[j - 1].position - target[i - 1].position,
                        *weight,
                    ));
                }
                index = (i - 1) * width + j - 1;
            }
            2 => index = (i - 1) * width + j,
            _ => index = i * width + j - 1,
        }
    }
    let (anchor, offset) = (!anchors.is_empty())
        .then(|| {
            let total_weight = anchors.iter().map(|(_, _, weight)| weight).sum::<f64>();
            let anchor = anchors
                .iter()
                .map(|(position, _, weight)| position * weight)
                .sum::<f64>()
                / total_weight;
            let offset = anchors
                .iter()
                .map(|(_, difference, weight)| difference * weight)
                .sum::<f64>()
                / total_weight;
            (anchor, offset)
        })
        .map_or((None, None), |(anchor, offset)| {
            (Some(anchor), Some(offset))
        });
    Alignment {
        score,
        anchor,
        offset,
    }
}

#[cfg(test)]
mod tests {
    use super::auto_locus_layouts;
    use crate::{Analysis, Cluster, Gene, GeneRef, Link, Locus};

    fn cluster(loci: usize, genes: usize) -> Cluster {
        Cluster {
            name: "test".into(),
            loci: (0..loci)
                .map(|locus| {
                    let start = locus * 100;
                    Locus {
                        name: "locus".into(),
                        start,
                        end: start + genes,
                        genes: (0..genes)
                            .map(|gene| Gene {
                                label: gene.to_string(),
                                names: Vec::new(),
                                start: start + gene,
                                end: start + gene + 1,
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
}
