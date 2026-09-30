use std::collections::BTreeMap;

use serde::Serialize;

use crate::layout::auto_locus_layouts;
use crate::{Analysis, GeneRef, Link, build_gene_groups};

/// JSON object consumed by clustermap.js.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotData {
    pub clusters: Vec<PlotCluster>,
    pub links: Vec<PlotLink>,
    pub groups: Vec<PlotGroup>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotCluster {
    pub uid: String,
    pub name: String,
    pub loci: Vec<PlotLocus>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotLocus {
    pub uid: String,
    pub name: String,
    pub start: usize,
    pub end: usize,
    pub genes: Vec<PlotGene>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotGene {
    pub uid: String,
    pub label: String,
    pub names: BTreeMap<String, String>,
    pub start: usize,
    pub end: usize,
    pub strand: i8,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotLink {
    pub uid: String,
    pub query: PlotGeneRef,
    pub target: PlotGeneRef,
    pub identity: f32,
    pub similarity: f32,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotGeneRef {
    pub uid: String,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotGroup {
    pub uid: String,
    pub label: String,
    pub genes: Vec<String>,
    pub hidden: bool,
    pub colour: Option<String>,
}

/// Convert analysis output to clustermap.js data contract
impl Analysis {
    pub fn to_plot_data(&self) -> PlotData {
        self.to_plot_data_in_order(&(0..self.clusters.len()).collect::<Vec<_>>())
    }

    /// Convert analysis output to plot data with clusters in `order`.
    /// Gene and link IDs retain their original analysis indices, so reordering
    /// clusters never invalidates cross-references in links or groups.
    pub fn to_plot_data_in_order(&self, order: &[usize]) -> PlotData {
        assert_eq!(order.len(), self.clusters.len());
        let clusters = order
            .iter()
            .copied()
            .map(|cluster_index| {
                let cluster = &self.clusters[cluster_index];
                PlotCluster {
                    uid: cluster_id(cluster_index),
                    name: cluster.name.clone(),
                    loci: cluster
                        .loci
                        .iter()
                        .enumerate()
                        .map(|(locus_index, locus)| PlotLocus {
                            uid: locus_id(cluster_index, locus_index),
                            name: locus.name.clone(),
                            start: locus.start,
                            end: locus.end,
                            genes: locus
                                .genes
                                .iter()
                                .enumerate()
                                .map(|(gene_index, gene)| PlotGene {
                                    uid: gene_id(cluster_index, locus_index, gene_index),
                                    label: gene.label.clone(),
                                    names: gene.names.iter().cloned().collect(),
                                    start: gene.start,
                                    end: gene.end,
                                    strand: gene.strand,
                                })
                                .collect(),
                        })
                        .collect(),
                }
            })
            .collect();

        let links = self
            .links
            .iter()
            .enumerate()
            .map(|(link_index, link)| PlotLink {
                uid: format!("link-{link_index}"),
                query: PlotGeneRef {
                    uid: gene_ref_id(link.query),
                },
                target: PlotGeneRef {
                    uid: gene_ref_id(link.target),
                },
                identity: link.identity,
                similarity: link.similarity,
            })
            .collect();

        PlotData {
            clusters,
            links,
            groups: build_plot_groups(&self.links),
        }
    }

    /// Emit plot data with loci reordered and flipped to follow homology anchors.
    pub fn to_auto_arranged_plot_data(&self, order: &[usize]) -> PlotData {
        let layouts = auto_locus_layouts(self, order);
        let mut data = self.to_plot_data_in_order(order);
        for (cluster, layout) in data.clusters.iter_mut().zip(layouts) {
            let original = std::mem::take(&mut cluster.loci);
            cluster.loci = layout
                .into_iter()
                .map(|placement| {
                    let mut locus = original[placement.locus].clone();
                    if placement.reversed {
                        let sum = locus.start + locus.end;
                        for gene in &mut locus.genes {
                            let (start, end) = (sum - gene.end, sum - gene.start);
                            gene.start = start;
                            gene.end = end;
                            gene.strand = -gene.strand;
                        }
                        locus.genes.reverse();
                    }
                    locus
                })
                .collect();
        }
        data
    }
}

/// Convert retained links into clustermap.js homology groups.
pub fn build_plot_groups(links: &[Link]) -> Vec<PlotGroup> {
    build_gene_groups(links)
        .into_iter()
        .enumerate()
        .map(|(index, group)| PlotGroup {
            uid: format!("group-{index}"),
            label: format!("Group {index}"),
            genes: group.genes.into_iter().map(gene_ref_id).collect(),
            hidden: false,
            colour: None,
        })
        .collect()
}

fn cluster_id(cluster: usize) -> String {
    format!("cluster-{cluster}")
}

fn locus_id(cluster: usize, locus: usize) -> String {
    format!("locus-{cluster}-{locus}")
}

fn gene_id(cluster: usize, locus: usize, gene: usize) -> String {
    format!("gene-{cluster}-{locus}-{gene}")
}

fn gene_ref_id(reference: GeneRef) -> String {
    gene_id(reference.cluster, reference.locus, reference.gene)
}

#[cfg(test)]
mod tests {
    use crate::{AnalysisOptions, InputFile, analyse_genbank};

    #[test]
    fn plot_data_uses_stable_ids_and_omits_translations() {
        let first =
            br#"LOCUS       FIRST                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag=\"first\"
ORIGIN
        1 atggcttaa
//
"#;
        let second =
            br#"LOCUS       SECOND                     9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag=\"second\"
ORIGIN
        1 atggcttaa
//
"#;
        let files = [
            InputFile {
                name: "first.gbk",
                bytes: first,
            },
            InputFile {
                name: "second.gbk",
                bytes: second,
            },
        ];

        let data = analyse_genbank(&files, AnalysisOptions::default())
            .unwrap()
            .to_plot_data();

        assert_eq!(data.clusters[0].uid, "cluster-0");
        assert_eq!(data.clusters[0].loci[0].genes[0].uid, "gene-0-0-0");
        assert_eq!(data.links[0].query.uid, "gene-0-0-0");
        assert_eq!(data.links[0].target.uid, "gene-1-0-0");
        assert_eq!(data.groups.len(), 1);
        assert_eq!(data.groups[0].genes, ["gene-0-0-0", "gene-1-0-0"]);
    }

    #[test]
    fn can_reorder_clusters_without_changing_link_ids() {
        let first =
            br#"LOCUS       FIRST                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag=\"first\"
ORIGIN
        1 atggcttaa
//
"#;
        let files = [
            InputFile {
                name: "first.gbk",
                bytes: first,
            },
            InputFile {
                name: "second.gbk",
                bytes: first,
            },
        ];
        let analysis = analyse_genbank(&files, AnalysisOptions::default()).unwrap();
        let data = analysis.to_plot_data_in_order(&[1, 0]);

        assert_eq!(data.clusters[0].uid, "cluster-1");
        assert_eq!(data.links[0].query.uid, "gene-0-0-0");
        assert_eq!(data.links[0].target.uid, "gene-1-0-0");
    }
}
