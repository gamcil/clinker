use clinker_core::{Analysis, AnalysisOptions, InputFile, analyse_genbank};

fn example_inputs() -> [InputFile<'static>; 5] {
    [
        InputFile {
            name: "P. vexata CBS 129021.gbk",
            bytes: include_bytes!("../../../../examples/P. vexata CBS 129021.gbk"),
        },
        InputFile {
            name: "A. versicolor CBS 583.65.gbk",
            bytes: include_bytes!("../../../../examples/A. versicolor CBS 583.65.gbk"),
        },
        InputFile {
            name: "A. mulundensis DSM 5745.gbk",
            bytes: include_bytes!("../../../../examples/A. mulundensis DSM 5745.gbk"),
        },
        InputFile {
            name: "A. burnettii MST-FP2249.gbk",
            bytes: include_bytes!("../../../../examples/A. burnettii MST-FP2249.gbk"),
        },
        InputFile {
            name: "A. alliaceus CBS 536.65.gbk",
            bytes: include_bytes!("../../../../examples/A. alliaceus CBS 536.65.gbk"),
        },
    ]
}

fn has_link(analysis: &Analysis, query: &str, target: &str) -> bool {
    analysis.links.iter().any(|link| {
        let query_label = &analysis.clusters[link.query.cluster].loci[link.query.locus].genes
            [link.query.gene]
            .label;
        let target_label = &analysis.clusters[link.target.cluster].loci[link.target.locus].genes
            [link.target.gene]
            .label;
        query_label == query && target_label == target
    })
}

/// Baseline generated with clinker 0.0.32 / Biopython 1.80 using its default
/// global BLOSUM62 configuration. This is intentionally ignored in ordinary
/// debug test runs: the exact all-protein comparison is a release-test/CI job.
#[test]
#[ignore = "runs all global protein comparisons across five real examples"]
fn matches_python_link_filtering_for_the_bundled_examples() {
    let files = example_inputs();

    let at_thirty_percent = analyse_genbank(
        &files,
        AnalysisOptions {
            identity_cutoff: 0.30,
        },
    )
    .unwrap();
    assert_eq!(at_thirty_percent.links.len(), 79);
    assert!(has_link(
        &at_thirty_percent,
        "OJI99004.1",
        "ncbi:ETB97_009431-T1"
    ));
    assert!(!has_link(
        &at_thirty_percent,
        "OJI99004.1",
        "ncbi:ETB97_008321-T1"
    ));
    assert!(
        at_thirty_percent
            .clusters
            .iter()
            .flat_map(|cluster| &cluster.loci)
            .flat_map(|locus| &locus.genes)
            .all(|gene| gene.start <= gene.end)
    );

    let at_seventy_percent = analyse_genbank(
        &files,
        AnalysisOptions {
            identity_cutoff: 0.70,
        },
    )
    .unwrap();
    assert_eq!(at_seventy_percent.links.len(), 37);
    assert!(has_link(
        &at_seventy_percent,
        "OJI99005.1",
        "XP_026607260.1"
    ));
    assert!(!has_link(
        &at_seventy_percent,
        "OJI99004.1",
        "ncbi:ETB97_009431-T1"
    ));
}
