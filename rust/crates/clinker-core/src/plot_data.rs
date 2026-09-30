use std::collections::BTreeMap;

use serde::Serialize;

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
        let clusters = self
            .clusters
            .iter()
            .enumerate()
            .map(|(cluster_index, cluster)| PlotCluster {
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
}
