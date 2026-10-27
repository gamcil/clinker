use std::collections::HashSet;

use clinker_core::{AnalysisOptions, InputFile, analyse_genbank};

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

    let analysis = analyse_genbank(&files, AnalysisOptions::default()).unwrap();

    assert_eq!(analysis.clusters.len(), files.len());
    assert!(analysis.links.iter().all(|link| {
        link.query.cluster != link.target.cluster
            && link.identity >= AnalysisOptions::default().identity_cutoff
            && link.similarity >= link.identity
    }));
    assert!(!analysis.links.is_empty());

    let order = analysis.cluster_order();
    assert_eq!(
        order.iter().copied().collect::<HashSet<_>>().len(),
        files.len()
    );

    let plot_data = analysis.to_auto_arranged_plot_data(&order);
    assert_eq!(plot_data.clusters.len(), files.len());
    assert_eq!(plot_data.links.len(), analysis.links.len());
    assert!(!plot_data.groups.is_empty());

    let json = serde_json::to_string(&plot_data).unwrap();
    assert!(!json.contains("translation"));
}
