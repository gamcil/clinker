use std::{
    fs,
    path::Path,
    process::{Command, Output},
};

use clinker_core::AnalysisSession;

fn run(directory: &Path, arguments: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_clinker"))
        .current_dir(directory)
        .args(arguments)
        .output()
        .unwrap()
}

fn fixture(directory: &Path) {
    fs::write(
        directory.join("input.gbk"),
        include_bytes!("fixtures/session.gbk"),
    )
    .unwrap();
}

fn create_session(directory: &Path) -> Vec<u8> {
    fixture(directory);
    let output = run(
        directory,
        &["input.gbk", "--no-align", "--session", "session.json"],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    fs::read(directory.join("session.json")).unwrap()
}

#[test]
fn invalid_plot_extension_leaves_the_session_and_other_outputs_unchanged() {
    let directory = tempfile::tempdir().unwrap();
    let previous = create_session(directory.path());
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--session",
            "session.json",
            "--output",
            "summary.txt",
            "--plot",
            "figure.invalid",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("--plot output must end"));
    assert_eq!(
        fs::read(directory.path().join("session.json")).unwrap(),
        previous
    );
    assert!(!directory.path().join("summary.txt").exists());
}

#[test]
fn an_existing_export_is_rejected_before_input_parsing_or_session_updates() {
    let directory = tempfile::tempdir().unwrap();
    let previous = create_session(directory.path());
    fs::write(directory.path().join("invalid.gbk"), "not GenBank").unwrap();
    fs::write(directory.path().join("figure.json"), "previous figure").unwrap();
    let output = run(
        directory.path(),
        &[
            "invalid.gbk",
            "--session",
            "session.json",
            "--plot",
            "figure.json",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("already exists"));
    assert_eq!(
        fs::read(directory.path().join("session.json")).unwrap(),
        previous
    );
    assert_eq!(
        fs::read_to_string(directory.path().join("figure.json")).unwrap(),
        "previous figure"
    );
}

#[test]
fn output_aliases_are_rejected_even_with_force() {
    let directory = tempfile::tempdir().unwrap();
    let previous = create_session(directory.path());
    fs::create_dir(directory.path().join("subdir")).unwrap();
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--session",
            "session.json",
            "--plot",
            "subdir/../session.json",
            "--force",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("more than one output"));
    assert_eq!(
        fs::read(directory.path().join("session.json")).unwrap(),
        previous
    );
}

#[test]
fn exports_cannot_overwrite_an_input_or_an_unmodified_session() {
    let directory = tempfile::tempdir().unwrap();
    let previous = create_session(directory.path());
    let input = fs::read(directory.path().join("input.gbk")).unwrap();
    for arguments in [
        vec!["input.gbk", "--output", "./input.gbk", "--force"],
        vec![
            "--session",
            "session.json",
            "--plot",
            "./session.json",
            "--force",
        ],
    ] {
        let output = run(directory.path(), &arguments);
        assert!(!output.status.success());
        assert!(String::from_utf8_lossy(&output.stderr).contains("also an input file"));
    }
    assert_eq!(
        fs::read(directory.path().join("session.json")).unwrap(),
        previous
    );
    assert_eq!(fs::read(directory.path().join("input.gbk")).unwrap(), input);
}

#[test]
fn invalid_output_directories_do_not_modify_the_session_or_publish_other_outputs() {
    let directory = tempfile::tempdir().unwrap();
    let previous = create_session(directory.path());
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--session",
            "session.json",
            "--plot",
            "figure.json",
            "--matrix-out",
            "missing/matrix.csv",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("output directory missing"));
    assert_eq!(
        fs::read(directory.path().join("session.json")).unwrap(),
        previous
    );
    assert!(!directory.path().join("figure.json").exists());
}

#[test]
fn force_replaces_exports_and_session_updates_need_no_force() {
    let directory = tempfile::tempdir().unwrap();
    create_session(directory.path());
    let files = ["summary.txt", "figure.json", "project.json", "matrix.csv"];
    for file in files {
        fs::write(directory.path().join(file), "old output").unwrap();
    }
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--session",
            "session.json",
            "--output",
            "summary.txt",
            "--plot",
            "figure.json",
            "--project",
            "project.json",
            "--matrix-out",
            "matrix.csv",
            "--force",
        ],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for file in files {
        assert_ne!(
            fs::read_to_string(directory.path().join(file)).unwrap(),
            "old output"
        );
    }
    assert!(
        serde_json::from_slice::<serde_json::Value>(
            &fs::read(directory.path().join("figure.json")).unwrap()
        )
        .is_ok()
    );
    assert!(
        clinker_core::ClinkerProject::from_json(
            &fs::read(directory.path().join("project.json")).unwrap()
        )
        .is_ok()
    );
    let output = run(
        directory.path(),
        &["input.gbk", "--no-align", "--session", "session.json"],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let session =
        AnalysisSession::from_json(&fs::read(directory.path().join("session.json")).unwrap())
            .unwrap();
    assert_eq!(session.analysis.clusters.len(), 3);
    assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 6);
}

#[test]
fn viewer_mode_rejects_output_flags_instead_of_ignoring_them() {
    let directory = tempfile::tempdir().unwrap();
    let output = run(
        directory.path(),
        &["--open-project", "saved.json", "--output", "summary.txt"],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("cannot be used with"));
}

#[test]
fn filesystem_case_aliases_do_not_bypass_output_or_input_protection() {
    let directory = tempfile::tempdir().unwrap();
    fixture(directory.path());
    if !directory.path().join("INPUT.GBK").exists() {
        return;
    }
    let previous = fs::read(directory.path().join("input.gbk")).unwrap();
    let output = run(
        directory.path(),
        &["input.gbk", "--output", "INPUT.GBK", "--force"],
    );
    assert!(!output.status.success());
    assert_eq!(
        fs::read(directory.path().join("input.gbk")).unwrap(),
        previous
    );
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--project",
            "Figure.json",
            "--plot",
            "figure.json",
            "--force",
        ],
    );
    assert!(!output.status.success());
    assert!(!directory.path().join("figure.json").exists());
    fs::create_dir(directory.path().join("exports")).unwrap();
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--project",
            "exports/Figure.json",
            "--plot",
            "EXPORTS/figure.json",
            "--force",
        ],
    );
    assert!(!output.status.success());
    assert!(!directory.path().join("exports/figure.json").exists());
}

#[test]
fn an_output_hard_link_to_an_input_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    fixture(directory.path());
    fs::hard_link(
        directory.path().join("input.gbk"),
        directory.path().join("summary.txt"),
    )
    .unwrap();
    let previous = fs::read(directory.path().join("input.gbk")).unwrap();
    let output = run(
        directory.path(),
        &["input.gbk", "--output", "summary.txt", "--force"],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("also an input file"));
    assert_eq!(
        fs::read(directory.path().join("summary.txt")).unwrap(),
        previous
    );
}

#[cfg(unix)]
#[test]
fn force_does_not_follow_or_replace_output_symlinks() {
    let directory = tempfile::tempdir().unwrap();
    fixture(directory.path());
    fs::write(directory.path().join("previous.txt"), "keep me").unwrap();
    std::os::unix::fs::symlink("previous.txt", directory.path().join("summary.txt")).unwrap();
    let output = run(
        directory.path(),
        &[
            "input.gbk",
            "--no-align",
            "--output",
            "summary.txt",
            "--force",
        ],
    );
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("must be a regular file"));
    assert!(
        fs::symlink_metadata(directory.path().join("summary.txt"))
            .unwrap()
            .is_symlink()
    );
    assert_eq!(
        fs::read_to_string(directory.path().join("previous.txt")).unwrap(),
        "keep me"
    );
}

#[cfg(unix)]
#[test]
fn an_unmodified_read_only_session_can_be_exported() {
    use std::os::unix::fs::PermissionsExt;
    let directory = tempfile::tempdir().unwrap();
    let previous = create_session(directory.path());
    fs::set_permissions(
        directory.path().join("session.json"),
        fs::Permissions::from_mode(0o400),
    )
    .unwrap();
    let output = run(
        directory.path(),
        &["--session", "session.json", "--plot", "figure.json"],
    );
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        fs::read(directory.path().join("session.json")).unwrap(),
        previous
    );
}
