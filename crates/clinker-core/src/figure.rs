use serde::Serialize;

use crate::plot_data::cluster_id;

use crate::{
    Analysis, ClusterOrdering, ClusterPairSimilarity, GeneGroupDefinition, PlotData,
    ValidationError,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RowOrder {
    Input,
    Similarity,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FigureStage {
    Similarity,
    Ordering,
    Groups,
    Layout,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct FigureProgress {
    pub stage: FigureStage,
    pub completed: usize,
    pub total: usize,
}

#[derive(Debug)]
pub struct PreparedFigure {
    pub plot_data: PlotData,
    pub similarity_matrix: Vec<Vec<ClusterPairSimilarity>>,
    pub ordering: ClusterOrdering,
    pub clustering_dendrogram: Option<ClusteringDendrogram>,
}

/// Renderer tree derived from the computed Ward hierarchy.
#[derive(Debug, Serialize)]
pub struct ClusteringDendrogram {
    root: String,
    nodes: Vec<ClusteringNode>,
    distance: DendrogramDistance,
}

#[derive(Debug, Serialize)]
#[serde(rename_all = "camelCase")]
struct ClusteringNode {
    uid: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    cluster: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    label: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    children: Option<[String; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    branch_length: Option<f64>,
}

#[derive(Debug, Serialize)]
struct DendrogramDistance {
    kind: &'static str,
    label: &'static str,
}

impl ClusteringDendrogram {
    fn from_ordering(analysis: &Analysis, ordering: &ClusterOrdering) -> Option<Self> {
        if ordering.merges.is_empty() {
            return None;
        }
        let count = analysis.clusters.len();
        let uid = |index| format!("clustering-{index}");
        let mut nodes: Vec<_> = analysis
            .clusters
            .iter()
            .enumerate()
            .map(|(index, cluster)| ClusteringNode {
                uid: uid(index),
                cluster: Some(cluster_id(index)),
                label: Some(cluster.name.clone()),
                children: None,
                branch_length: None,
            })
            .collect();
        for (index, children) in ordering.merges.iter().enumerate() {
            let height = ordering.heights[index];
            for &child in children {
                let child_height = if child < count {
                    0.0
                } else {
                    ordering.heights[child - count]
                };
                nodes[child].branch_length = Some((height - child_height).max(0.0));
            }
            nodes.push(ClusteringNode {
                uid: uid(count + index),
                cluster: None,
                label: None,
                children: Some(children.map(uid)),
                branch_length: None,
            });
        }
        Some(Self {
            root: uid(nodes.len() - 1),
            nodes,
            distance: DendrogramDistance {
                kind: "mergeHeight",
                label: "Ward merge height",
            },
        })
    }
}

impl Analysis {
    /// Score, order, group and arrange a figure without changing the analysis.
    /// Similarity rows and hierarchy leaves refer to original cluster indices.
    pub fn prepare_figure(
        &self,
        row_order: RowOrder,
        definitions: &[GeneGroupDefinition],
        mut progress: impl FnMut(FigureProgress),
    ) -> Result<PreparedFigure, ValidationError> {
        let count = self.clusters.len();
        let total_pairs = count.saturating_mul(count.saturating_sub(1)) / 2;
        let similarity_matrix = self.cluster_similarity_matrix_with_progress(|completed| {
            progress(FigureProgress {
                stage: FigureStage::Similarity,
                completed,
                total: total_pairs,
            });
        })?;
        progress(FigureProgress {
            stage: FigureStage::Ordering,
            completed: 0,
            total: count,
        });
        let ordering = match row_order {
            RowOrder::Input => ClusterOrdering {
                order: (0..count).collect(),
                merges: Vec::new(),
                heights: Vec::new(),
            },
            RowOrder::Similarity => {
                self.cluster_ordering_from_similarity_matrix(&similarity_matrix)
            }
        };
        progress(FigureProgress {
            stage: FigureStage::Ordering,
            completed: count,
            total: count,
        });
        let plot_data = self.arrange_plot_data(&ordering.order, definitions, &mut progress);
        let clustering_dendrogram = ClusteringDendrogram::from_ordering(self, &ordering);
        Ok(PreparedFigure {
            plot_data,
            similarity_matrix,
            ordering,
            clustering_dendrogram,
        })
    }
}
