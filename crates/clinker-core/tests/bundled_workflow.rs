use std::collections::HashSet;

use clinker_core::{AnalysisOptions, FigureStage, InputFile, RowOrder, analyse};

#[test]
fn bundled_genbank_files_complete_the_analysis_and_plot_data_workflow() {
    let files = [
        InputFile {
            name: "A. versicolor CBS 583.65.gbk",
            bytes: include_bytes!("fixtures/examples/A. versicolor CBS 583.65.gbk"),
        },
        InputFile {
            name: "A. burnettii MST-FP2249.gbk",
            bytes: include_bytes!("fixtures/examples/A. burnettii MST-FP2249.gbk"),
        },
    ];

    let analysis = analyse(&files, AnalysisOptions::default()).unwrap();

    assert_eq!(analysis.clusters.len(), files.len());
    assert!(analysis.links.iter().all(|link| {
        link.query.cluster != link.target.cluster
            && link.identity >= AnalysisOptions::default().identity_cutoff
            && link.similarity >= link.identity
    }));
    assert!(!analysis.links.is_empty());

    let mut updates = Vec::new();
    let figure = analysis
        .prepare_figure(RowOrder::Similarity, &[], |update| updates.push(update))
        .unwrap();
    let order = &figure.ordering.order;
    assert_eq!(figure.ordering.merges.len(), files.len() - 1);
    for stage in [
        FigureStage::Similarity,
        FigureStage::Ordering,
        FigureStage::Groups,
        FigureStage::Layout,
    ] {
        let last = updates
            .iter()
            .rfind(|update| update.stage == stage)
            .unwrap();
        assert_eq!(last.completed, last.total, "{stage:?}");
    }
    assert!(
        updates
            .iter()
            .all(|update| update.completed <= update.total)
    );
    assert_eq!(
        order.iter().copied().collect::<HashSet<_>>().len(),
        files.len()
    );

    let plot_data = figure.plot_data;
    let input = analysis
        .prepare_figure(RowOrder::Input, &[], |_| {})
        .unwrap();
    assert_eq!(input.ordering.order, [0, 1]);
    assert!(input.ordering.merges.is_empty());
    assert!(input.ordering.heights.is_empty());
    assert_eq!(input.similarity_matrix, figure.similarity_matrix);
    assert_eq!(plot_data.clusters.len(), files.len());
    assert_eq!(plot_data.links.len(), analysis.links.len());
    assert!(!plot_data.groups.is_empty());

    let json = serde_json::to_string(&plot_data).unwrap();
    assert!(!json.contains("translation"));
}
