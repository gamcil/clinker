use std::{
    fs,
    path::Path,
    process::{Command, Output},
};

use clinker_core::{AlignmentSettings, AnalysisSession, ClinkerProject, KmerPrefilter};

fn run(directory: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_clinker"))
        .current_dir(directory)
        .args(args)
        .output()
        .unwrap()
}

fn succeeds(directory: &Path, args: &[&str]) -> Output {
    let output = run(directory, args);
    assert!(
        output.status.success(),
        "{args:?}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output
}

fn session(directory: &Path) -> AnalysisSession {
    AnalysisSession::from_json(&fs::read(directory.join("session.json")).unwrap()).unwrap()
}

fn genbank(directory: &Path) {
    fs::write(
        directory.join("input.gbk"),
        include_bytes!("fixtures/session.gbk"),
    )
    .unwrap();
}

fn gff(directory: &Path, sequence: &str) {
    fs::create_dir_all(directory).unwrap();
    fs::write(
        directory.join("sample.gff3"),
        "##gff-version 3\nseq\t.\tCDS\t1\t9\t.\t+\t0\tID=gene\n",
    )
    .unwrap();
    fs::write(directory.join("sample.fna"), format!(">seq\n{sequence}\n")).unwrap();
}

#[test]
fn gff_discovery_preserves_directory_pairing_order_and_distinct_inputs() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    gff(&root.join("a"), "ATGGCTTAA");
    gff(&root.join("b"), "ATGAAATAA");
    succeeds(
        root,
        &[
            "b",
            "a",
            "./a/sample.gff3",
            "--no-align",
            "--session",
            "session.json",
        ],
    );
    let analysis = session(root).analysis;
    assert_eq!(analysis.clusters.len(), 2);
    assert_eq!(analysis.clusters[0].loci[0].genes[0].translation, "MK");
    assert_eq!(analysis.clusters[1].loci[0].genes[0].translation, "MA");
    succeeds(
        &root.join("a"),
        &["sample.gff3", "--no-align", "--plot", "figure.json"],
    );
    let figure: serde_json::Value =
        serde_json::from_slice(&fs::read(root.join("a/figure.json")).unwrap()).unwrap();
    assert_eq!(figure["clusters"][0]["source_file"], "sample.gff3");
}

#[test]
fn ambiguous_or_orphan_fastas_fail_and_explicit_selection_resolves_the_pair() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    gff(root, "ATGGCTTAA");
    fs::write(root.join("sample.fa"), ">seq\nATGAAATAA\n").unwrap();
    fs::write(root.join("orphan.fna"), ">seq\nATGAAATAA\n").unwrap();
    for args in [
        vec!["sample.gff3", "--no-align"],
        vec!["sample.gff3", "orphan.fna", "--no-align"],
    ] {
        assert!(!run(root, &args).status.success(), "{args:?}");
    }
    succeeds(
        root,
        &[
            "sample.gff3",
            "sample.fa",
            "--no-align",
            "--session",
            "session.json",
        ],
    );
    assert_eq!(
        session(root).analysis.clusters[0].loci[0].genes[0].translation,
        "MK"
    );
}

#[test]
fn sessions_extend_and_export_restorable_projects_figures_and_escaped_csv() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    genbank(root);
    succeeds(
        root,
        &["input.gbk", "--no-align", "--session", "session.json"],
    );
    succeeds(
        root,
        &["input.gbk", "--session", "session.json", "--jobs", "2"],
    );
    assert_eq!(session(root).alignment(), Some(AlignmentSettings::Skipped));
    assert!(session(root).analysis.links.is_empty());
    succeeds(root, &["--session", "session.json", "--reanalyse"]);
    let mut saved = session(root);
    assert_eq!(saved.analysis.clusters.len(), 2);
    assert!(!saved.analysis.links.is_empty());
    let aligned = succeeds(root, &["--session", "session.json", "--decimals", "4"]);
    assert_eq!(
        String::from_utf8(aligned.stdout).unwrap(),
        "input vs input\n--------------\nQuery  Target  Identity  Similarity\ngene   gene    1.0000    1.0000    \n"
    );
    let delimited = succeeds(
        root,
        &[
            "--session",
            "session.json",
            "--delimiter",
            ",",
            "--hide-alignment-headers",
            "--hide-link-headers",
        ],
    );
    assert_eq!(
        String::from_utf8(delimited.stdout).unwrap(),
        "gene,gene,1.00,1.00\n"
    );
    saved.analysis.clusters[0].name = "one, with comma".into();
    saved.analysis.clusters[1].name = "two \"quoted\"\nsecond line".into();
    fs::write(root.join("session.json"), saved.to_json().unwrap()).unwrap();
    succeeds(
        root,
        &[
            "--session",
            "session.json",
            "--output",
            "summary.txt",
            "--plot",
            "figure.html",
            "--project",
            "project.json",
            "--matrix-out",
            "matrix.csv",
        ],
    );
    let project = ClinkerProject::from_json(&fs::read(root.join("project.json")).unwrap()).unwrap();
    assert_eq!(project.analysis, saved);
    assert_eq!(project.data["clusters"].as_array().unwrap().len(), 2);
    assert!(
        fs::read_to_string(root.join("figure.html"))
            .unwrap()
            .starts_with("<!doctype html>")
    );
    assert!(
        fs::read_to_string(root.join("summary.txt"))
            .unwrap()
            .contains("Identity")
    );
    let mut csv = csv::Reader::from_path(root.join("matrix.csv")).unwrap();
    assert_eq!(csv.headers().unwrap().get(1), Some("one, with comma"));
    assert_eq!(
        csv.records().nth(1).unwrap().unwrap().get(0),
        Some("two \"quoted\"\nsecond line")
    );
    succeeds(
        root,
        &[
            "--session",
            "session.json",
            "--plot",
            "figure.html",
            "--force",
        ],
    );
    assert_eq!(session(root), saved);
}

#[test]
fn session_settings_are_inherited_and_changes_rebuild_all_comparisons() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    genbank(root);
    succeeds(
        root,
        &[
            "input.gbk",
            "--session",
            "session.json",
            "--prefilter",
            "--kmer-size",
            "5",
            "--min-shared-kmers",
            "8",
        ],
    );
    succeeds(root, &["input.gbk", "--session", "session.json"]);
    let saved = session(root);
    assert_eq!(
        saved.alignment(),
        Some(AlignmentSettings::Aligned {
            prefilter: Some(KmerPrefilter {
                kmer_size: 5,
                min_shared_kmers: 8
            }),
        })
    );
    // This low-complexity protein has only two distinct words: filtering removes its match.
    assert!(saved.analysis.links.is_empty());
    let previous = fs::read(root.join("session.json")).unwrap();
    for extra in [
        &["--identity", "0.8"][..],
        &["--no-prefilter"][..],
        &["--no-align"][..],
    ] {
        let mut args = vec!["input.gbk", "--session", "session.json"];
        args.extend(extra);
        let result = run(root, &args);
        assert!(!result.status.success());
        assert!(String::from_utf8_lossy(&result.stderr).contains("--reanalyse"));
        assert_eq!(fs::read(root.join("session.json")).unwrap(), previous);
    }
    succeeds(
        root,
        &[
            "--session",
            "session.json",
            "--reanalyse",
            "--no-prefilter",
            "--identity",
            "0.8",
        ],
    );
    let rebuilt = session(root);
    assert_eq!(rebuilt.analysis.clusters.len(), 2);
    assert_eq!(rebuilt.analysis.links.len(), 1);
    assert_eq!(rebuilt.identity_cutoff(), 0.8);
    assert_eq!(
        rebuilt.alignment(),
        Some(AlignmentSettings::Aligned { prefilter: None })
    );
    succeeds(
        root,
        &["input.gbk", "--session", "session.json", "--jobs", "2"],
    );
    assert_eq!(session(root).analysis.links.len(), 3);

    let mut legacy = serde_json::to_value(session(root)).unwrap();
    legacy.as_object_mut().unwrap().remove("alignment");
    fs::write(
        root.join("session.json"),
        serde_json::to_vec(&legacy).unwrap(),
    )
    .unwrap();
    succeeds(
        root,
        &["--session", "session.json", "--project", "legacy.json"],
    );
    assert_eq!(session(root).alignment(), None);
    let rejected = run(root, &["input.gbk", "--session", "session.json"]);
    assert!(!rejected.status.success());
    assert!(String::from_utf8_lossy(&rejected.stderr).contains("--reanalyse"));
    succeeds(root, &["--session", "session.json", "--reanalyse"]);
    assert_eq!(session(root).analysis.links.len(), 3);
    assert!(session(root).alignment().is_some());
    succeeds(
        root,
        &["--session", "session.json", "--reanalyse", "--no-align"],
    );
    assert_eq!(session(root).alignment(), Some(AlignmentSettings::Skipped));
    assert!(session(root).analysis.links.is_empty());
}

#[test]
fn rejected_and_failed_runs_preserve_inputs_sessions_and_previous_exports() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    genbank(root);
    succeeds(
        root,
        &["input.gbk", "--no-align", "--session", "session.json"],
    );
    fs::write(root.join("bad.gbk"), "invalid GenBank").unwrap();
    fs::write(root.join("figure.json"), "previous figure").unwrap();
    fs::create_dir(root.join("subdir")).unwrap();
    let previous_session = fs::read(root.join("session.json")).unwrap();
    let previous_input = fs::read(root.join("input.gbk")).unwrap();
    let cases: &[&[&str]] = &[
        &[
            "input.gbk",
            "--session",
            "session.json",
            "--plot",
            "figure.invalid",
        ],
        &[
            "bad.gbk",
            "--session",
            "session.json",
            "--plot",
            "figure.json",
        ],
        &[
            "bad.gbk",
            "--session",
            "session.json",
            "--plot",
            "figure.json",
            "--force",
        ],
        &[
            "input.gbk",
            "--session",
            "session.json",
            "--plot",
            "subdir/../session.json",
            "--force",
        ],
        &[
            "--session",
            "session.json",
            "--plot",
            "session.json",
            "--force",
        ],
        &["input.gbk", "--output", "./input.gbk", "--force"],
        &[
            "input.gbk",
            "--session",
            "session.json",
            "--matrix-out",
            "missing/matrix.csv",
        ],
        &["--open-project", "project.json", "--output", "summary.txt"],
    ];
    for args in cases {
        let result = run(root, args);
        assert!(!result.status.success(), "{args:?}");
        assert_eq!(
            fs::read(root.join("session.json")).unwrap(),
            previous_session,
            "{args:?}"
        );
        assert_eq!(
            fs::read(root.join("input.gbk")).unwrap(),
            previous_input,
            "{args:?}"
        );
        assert_eq!(
            fs::read_to_string(root.join("figure.json")).unwrap(),
            "previous figure",
            "{args:?}"
        );
    }
    assert_eq!(
        fs::read_dir(root).unwrap().count(),
        5,
        "temporary outputs are cleaned up"
    );
}

#[test]
fn filesystem_aliases_cannot_bypass_input_or_output_protection() {
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    genbank(root);
    let original = fs::read(root.join("input.gbk")).unwrap();
    fs::hard_link(root.join("input.gbk"), root.join("summary.txt")).unwrap();
    assert!(
        !run(root, &["input.gbk", "--output", "summary.txt", "--force"])
            .status
            .success()
    );
    if root.join("INPUT.GBK").exists() {
        for args in [
            vec!["input.gbk", "--output", "INPUT.GBK", "--force"],
            vec![
                "input.gbk",
                "--project",
                "Figure.json",
                "--plot",
                "figure.json",
                "--force",
            ],
        ] {
            assert!(!run(root, &args).status.success(), "{args:?}");
        }
        assert!(!root.join("figure.json").exists());
    }
    assert_eq!(fs::read(root.join("input.gbk")).unwrap(), original);
}

#[cfg(unix)]
#[test]
fn symlinks_do_not_break_pairing_recurse_forever_or_redirect_exports() {
    use std::os::unix::fs::symlink;
    let directory = tempfile::tempdir().unwrap();
    let root = directory.path();
    gff(&root.join("inputs"), "ATGGCTTAA");
    fs::rename(
        root.join("inputs/sample.fna"),
        root.join("inputs/sequence.txt"),
    )
    .unwrap();
    symlink("sequence.txt", root.join("inputs/sample.fna")).unwrap();
    symlink(".", root.join("inputs/loop")).unwrap();
    succeeds(root, &["inputs", "--no-align", "--session", "session.json"]);
    assert_eq!(
        session(root).analysis.clusters[0].loci[0].genes[0].translation,
        "MA"
    );
    symlink("session.json", root.join("figure.json")).unwrap();
    assert!(
        !run(
            root,
            &[
                "--session",
                "session.json",
                "--plot",
                "figure.json",
                "--force"
            ]
        )
        .status
        .success()
    );
    assert!(
        fs::symlink_metadata(root.join("figure.json"))
            .unwrap()
            .is_symlink()
    );
}
