use std::{
    fs,
    path::{Path, PathBuf},
    process::Command,
    time::{SystemTime, UNIX_EPOCH},
};

use clinker_core::AnalysisSession;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../../examples")
        .join(name)
}

fn run(arguments: &[&Path]) {
    let output = Command::new(env!("CARGO_BIN_EXE_clinker-rs"))
        .args(arguments)
        .output()
        .expect("run clinker CLI");
    assert!(
        output.status.success(),
        "CLI failed:\nstdout:\n{}\nstderr:\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr),
    );
}

#[test]
fn creates_reopens_and_extends_an_analysis_session() {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .unwrap()
        .as_nanos();
    let directory = std::env::temp_dir().join(format!(
        "clinker-session-test-{}-{unique}",
        std::process::id()
    ));
    fs::create_dir(&directory).unwrap();
    let session_path = directory.join("analysis.json");
    let plot_path = directory.join("plot.json");
    let first = fixture("A. versicolor CBS 583.65.gbk");
    let second = fixture("A. alliaceus CBS 536.65.gbk");

    run(&[
        first.as_path(),
        Path::new("--no-align"),
        Path::new("--session"),
        session_path.as_path(),
    ]);
    let initial = AnalysisSession::from_json(&fs::read(&session_path).unwrap()).unwrap();
    assert_eq!(initial.analysis.clusters.len(), 1);
    assert!(initial.analysis.links.is_empty());

    run(&[
        second.as_path(),
        Path::new("--session"),
        session_path.as_path(),
    ]);
    let extended = AnalysisSession::from_json(&fs::read(&session_path).unwrap()).unwrap();
    assert_eq!(extended.analysis.clusters.len(), 2);
    assert!(!extended.analysis.links.is_empty());

    run(&[
        Path::new("--session"),
        session_path.as_path(),
        Path::new("--plot"),
        plot_path.as_path(),
    ]);
    assert!(plot_path.exists());

    fs::remove_dir_all(directory).unwrap();
}
