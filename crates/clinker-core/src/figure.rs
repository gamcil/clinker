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
        self.validate()?;
        let count = self.clusters.len();
        let total_pairs = count.saturating_mul(count.saturating_sub(1)) / 2;
        let similarity_matrix = self.cluster_similarity_matrix_with_progress(|completed| {
            progress(FigureProgress {
                stage: FigureStage::Similarity,
                completed,
                total: total_pairs,
            });
        });
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
        Ok(PreparedFigure {
            plot_data,
            similarity_matrix,
            ordering,
        })
    }
}
