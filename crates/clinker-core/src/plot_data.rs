use std::collections::{BTreeMap, BTreeSet};

use serde::{Deserialize, Serialize};

use crate::layout::auto_locus_layouts_with_progress;
use crate::{
    Analysis, FigureProgress, FigureStage, GeneRef, ValidationError, build_gene_groups,
    build_gene_groups_with_progress,
};

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
    #[serde(skip_serializing_if = "Option::is_none")]
    pub source_file: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub organism: Option<String>,
    pub loci: Vec<PlotLocus>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotLocus {
    pub uid: String,
    pub name: String,
    pub start: u64,
    pub end: u64,
    /// Initial horizontal placement in sequence-coordinate units.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub offset: Option<f64>,
    pub genes: Vec<PlotGene>,
}

#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct PlotGene {
    pub uid: String,
    pub label: String,
    pub names: BTreeMap<String, String>,
    pub start: u64,
    pub end: u64,
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

/// A user-defined functional group. Gene names are matched against both the
/// displayed gene label and retained qualifier values.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize, Serialize)]
pub struct GeneGroupDefinition {
    pub label: String,
    pub genes: Vec<String>,
    pub colour: Option<String>,
}

/// Convert analysis output to clustermap.js data contract
impl Analysis {
    pub fn to_plot_data(&self) -> Result<PlotData, ValidationError> {
        self.to_plot_data_in_order(&(0..self.clusters.len()).collect::<Vec<_>>())
    }

    pub fn to_plot_data_with_group_definitions(
        &self,
        definitions: &[GeneGroupDefinition],
    ) -> Result<PlotData, ValidationError> {
        self.to_plot_data_in_order_with_group_definitions(
            &(0..self.clusters.len()).collect::<Vec<_>>(),
            definitions,
        )
    }

    /// Convert analysis output to plot data with clusters in `order`.
    /// Gene and link IDs retain their original analysis indices, so reordering
    /// clusters never invalidates cross-references in links or groups.
    pub fn to_plot_data_in_order(&self, order: &[usize]) -> Result<PlotData, ValidationError> {
        self.to_plot_data_in_order_with_group_definitions(order, &[])
    }

    pub fn to_plot_data_in_order_with_group_definitions(
        &self,
        order: &[usize],
        definitions: &[GeneGroupDefinition],
    ) -> Result<PlotData, ValidationError> {
        self.validate()?;
        self.validate_order(order)?;
        Ok(self.to_plot_data_in_order_with_groups(order, build_plot_groups(self, definitions)))
    }

    fn to_plot_data_in_order_with_groups(
        &self,
        order: &[usize],
        groups: Vec<PlotGroup>,
    ) -> PlotData {
        let clusters = order
            .iter()
            .copied()
            .map(|cluster_index| {
                let cluster = &self.clusters[cluster_index];
                PlotCluster {
                    uid: cluster_id(cluster_index),
                    name: cluster.name.clone(),
                    source_file: cluster.source_file.clone(),
                    organism: cluster.organism.clone(),
                    loci: cluster
                        .loci
                        .iter()
                        .enumerate()
                        .map(|(locus_index, locus)| PlotLocus {
                            uid: locus_id(cluster_index, locus_index),
                            name: locus.name.clone(),
                            start: locus.start,
                            end: locus.end,
                            offset: None,
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
            groups,
        }
    }

    /// Emit plot data with loci reordered and flipped to follow homology anchors.
    pub fn to_auto_arranged_plot_data(&self, order: &[usize]) -> Result<PlotData, ValidationError> {
        self.to_auto_arranged_plot_data_with_group_definitions(order, &[])
    }

    pub fn to_auto_arranged_plot_data_with_group_definitions(
        &self,
        order: &[usize],
        definitions: &[GeneGroupDefinition],
    ) -> Result<PlotData, ValidationError> {
        self.validate()?;
        self.validate_order(order)?;
        Ok(self.arrange_plot_data(order, definitions, &mut |_| {}))
    }

    pub(crate) fn arrange_plot_data(
        &self,
        order: &[usize],
        definitions: &[GeneGroupDefinition],
        progress: &mut impl FnMut(FigureProgress),
    ) -> PlotData {
        let groups = build_plot_groups_with_progress(self, definitions, |completed| {
            progress(FigureProgress {
                stage: FigureStage::Groups,
                completed,
                total: self.links.len().saturating_mul(2),
            });
        });
        let layouts = auto_locus_layouts_with_progress(self, order, |completed| {
            progress(FigureProgress {
                stage: FigureStage::Layout,
                completed,
                total: self.clusters.len(),
            });
        });
        let mut data = self.to_plot_data_in_order_with_groups(order, groups);
        for (cluster, layout) in data.clusters.iter_mut().zip(layouts) {
            let mut original = std::mem::take(&mut cluster.loci)
                .into_iter()
                .map(Some)
                .collect::<Vec<_>>();
            cluster.loci = layout
                .into_iter()
                .map(|placement| {
                    let mut locus = original[placement.locus]
                        .take()
                        .expect("layout contains each locus exactly once");
                    locus.offset = placement.offset;
                    if placement.reversed {
                        for gene in &mut locus.genes {
                            let (start, end) = (
                                locus.start + (locus.end - gene.end),
                                locus.start + (locus.end - gene.start),
                            );
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

/// Convert configured functions and retained links into clustermap.js groups.
pub fn build_plot_groups(
    analysis: &Analysis,
    definitions: &[GeneGroupDefinition],
) -> Vec<PlotGroup> {
    configured_plot_groups(analysis, definitions, build_gene_groups(&analysis.links))
}

fn build_plot_groups_with_progress(
    analysis: &Analysis,
    definitions: &[GeneGroupDefinition],
    progress: impl FnMut(usize),
) -> Vec<PlotGroup> {
    configured_plot_groups(
        analysis,
        definitions,
        build_gene_groups_with_progress(&analysis.links, progress),
    )
}

fn configured_plot_groups(
    analysis: &Analysis,
    definitions: &[GeneGroupDefinition],
    automatic: Vec<crate::GeneGroup>,
) -> Vec<PlotGroup> {
    let mut groups = definitions
        .iter()
        .map(|definition| ConfiguredGroup {
            label: definition.label.clone(),
            genes: BTreeSet::new(),
            colour: definition.colour.clone(),
        })
        .collect::<Vec<_>>();

    let mut requested = BTreeMap::<&str, Vec<usize>>::new();
    for (index, definition) in definitions.iter().enumerate() {
        for name in &definition.genes {
            requested.entry(name).or_default().push(index);
        }
    }
    if !requested.is_empty() {
        for (cluster, cluster_data) in analysis.clusters.iter().enumerate() {
            for (locus, locus_data) in cluster_data.loci.iter().enumerate() {
                for (gene, gene_data) in locus_data.genes.iter().enumerate() {
                    for name in std::iter::once(gene_data.label.as_str())
                        .chain(gene_data.names.iter().map(|(_, value)| value.as_str()))
                    {
                        if let Some(indices) = requested.get(name) {
                            for &index in indices {
                                groups[index].genes.insert(GeneRef {
                                    cluster,
                                    locus,
                                    gene,
                                });
                            }
                        }
                    }
                }
            }
        }
    }

    for automatic_group in automatic {
        let genes = automatic_group.genes.into_iter().collect::<BTreeSet<_>>();
        let overlaps = groups
            .iter()
            .enumerate()
            .filter_map(|(index, group)| (!group.genes.is_disjoint(&genes)).then_some(index))
            .collect::<Vec<_>>();
        if let Some((&keep, merged)) = overlaps.split_first() {
            groups[keep].genes.extend(genes);
            for &index in merged.iter().rev() {
                let removed = groups.remove(index);
                groups[keep].genes.extend(removed.genes);
            }
        } else {
            groups.push(ConfiguredGroup {
                label: String::new(),
                genes,
                colour: None,
            });
        }
    }

    groups
        .into_iter()
        .filter(|group| !group.genes.is_empty())
        .enumerate()
        .map(|(index, group)| PlotGroup {
            uid: format!("group-{index}"),
            label: if group.label.is_empty() {
                format!("Group {index}")
            } else {
                group.label
            },
            genes: group.genes.into_iter().map(gene_ref_id).collect(),
            hidden: false,
            colour: group.colour,
        })
        .collect()
}

struct ConfiguredGroup {
    label: String,
    genes: BTreeSet<GeneRef>,
    colour: Option<String>,
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
    use crate::{AnalysisOptions, GeneGroupDefinition, InputFile, analyse_genbank};

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
            .to_plot_data()
            .unwrap();

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
        let data = analysis.to_plot_data_in_order(&[1, 0]).unwrap();

        for order in [&[0, 0][..], &[0, 2][..], &[0][..]] {
            assert!(analysis.to_plot_data_in_order(order).is_err());
            assert!(analysis.to_auto_arranged_plot_data(order).is_err());
        }

        assert_eq!(data.clusters[0].uid, "cluster-1");
        assert_eq!(data.links[0].query.uid, "gene-0-0-0");
        assert_eq!(data.links[0].target.uid, "gene-1-0-0");
    }

    #[test]
    fn auto_layout_emits_locus_offsets_in_plot_data() {
        let input =
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
                bytes: input,
            },
            InputFile {
                name: "second.gbk",
                bytes: input,
            },
        ];
        let analysis = analyse_genbank(&files, AnalysisOptions::default()).unwrap();

        let data = analysis.to_auto_arranged_plot_data(&[0, 1]).unwrap();

        // A zero offset is meaningful: it overrides default locus packing and
        // aligns this homologous locus to the first cluster.
        assert_eq!(data.clusters[1].loci[0].offset, Some(0.0));
    }

    #[test]
    fn function_groups_seed_link_groups_and_work_without_alignment() {
        let first =
            br#"LOCUS       FIRST                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="first"
ORIGIN
        1 atggcttaa
//
"#;
        let second =
            br#"LOCUS       SECOND                     9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="second"
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
        let mut analysis = analyse_genbank(&files, AnalysisOptions::default()).unwrap();
        let definitions = [
            GeneGroupDefinition {
                label: "Core enzyme".into(),
                genes: vec!["first".into()],
                colour: Some("#336699".into()),
            },
            GeneGroupDefinition {
                label: "Other annotation".into(),
                genes: vec!["second".into()],
                colour: Some("#993366".into()),
            },
        ];

        let linked = analysis
            .to_plot_data_with_group_definitions(&definitions)
            .unwrap();
        assert_eq!(linked.groups.len(), 1);
        assert_eq!(linked.groups[0].label, "Core enzyme");
        assert_eq!(linked.groups[0].colour.as_deref(), Some("#336699"));
        assert_eq!(linked.groups[0].genes.len(), 2);

        analysis.links.clear();
        let annotation_only = analysis
            .to_plot_data_with_group_definitions(&definitions)
            .unwrap();
        assert!(annotation_only.links.is_empty());
        assert_eq!(annotation_only.groups.len(), 2);
        assert_eq!(annotation_only.groups[0].genes, ["gene-0-0-0"]);
        assert_eq!(annotation_only.groups[1].genes, ["gene-1-0-0"]);
    }
}
