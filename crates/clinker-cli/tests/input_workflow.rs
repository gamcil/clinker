use std::{fs, path::Path, process::Command};

use clinker_core::AnalysisSession;

const GFF: &str = "##gff-version 3\nseq\t.\tCDS\t1\t9\t.\t+\t0\tID=gene\n";

fn write_pair(directory: &Path, fasta: &str) {
    fs::create_dir_all(directory).unwrap();
    fs::write(directory.join("sample.gff3"), GFF).unwrap();
    fs::write(directory.join("sample.fna"), fasta).unwrap();
}

fn analyse(directory: &Path, inputs: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_clinker"))
        .current_dir(directory)
        .args(inputs)
        .args(["--no-align", "--session", "analysis.json"])
        .output()
        .unwrap()
}

fn session(directory: &Path) -> AnalysisSession {
    AnalysisSession::from_json(&fs::read(directory.join("analysis.json")).unwrap()).unwrap()
}

#[test]
fn bare_gff_filename_discovers_its_companion() {
    let directory = tempfile::tempdir().unwrap();
    write_pair(directory.path(), ">seq\nATGGCTTAA\n");
    let output = analyse(directory.path(), &["sample.gff3"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        session(directory.path()).analysis.clusters[0].loci[0].genes[0].translation,
        "MA"
    );
}

#[test]
fn duplicate_basenames_are_paired_within_their_own_directories_in_input_order() {
    let directory = tempfile::tempdir().unwrap();
    write_pair(&directory.path().join("a"), ">seq\nATGGCTTAA\n");
    write_pair(&directory.path().join("b"), ">seq\nATGAAATAA\n");
    let output = analyse(directory.path(), &["b", "a"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let session = session(directory.path());
    let clusters = &session.analysis.clusters;
    assert_eq!(clusters.len(), 2);
    assert_eq!(clusters[0].source_file.as_deref(), Some("sample.gff3"));
    assert_eq!(clusters[1].source_file.as_deref(), Some("sample.gff3"));
    assert_eq!(clusters[0].loci[0].genes[0].translation, "MK");
    assert_eq!(clusters[1].loci[0].genes[0].translation, "MA");
}

#[test]
fn ambiguous_companions_require_an_explicit_selection() {
    let directory = tempfile::tempdir().unwrap();
    write_pair(directory.path(), ">seq\nATGGCTTAA\n");
    fs::write(directory.path().join("sample.fa"), ">seq\nATGAAATAA\n").unwrap();
    let output = analyse(directory.path(), &["sample.gff3"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("ambiguous FASTA companions"));
    assert!(!directory.path().join("analysis.json").exists());
    let output = analyse(directory.path(), &["sample.gff3", "sample.fa"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        session(directory.path()).analysis.clusters[0].loci[0].genes[0].translation,
        "MK"
    );
}

#[test]
fn unmatched_explicit_fasta_is_rejected() {
    let directory = tempfile::tempdir().unwrap();
    write_pair(directory.path(), ">seq\nATGGCTTAA\n");
    fs::write(directory.path().join("other.fna"), ">seq\nATGAAATAA\n").unwrap();
    let output = analyse(directory.path(), &["sample.gff3", "other.fna"]);
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("no GFF3/GFF input"));
}

#[test]
fn discovery_sorts_directory_entries_and_deduplicates_path_aliases() {
    let directory = tempfile::tempdir().unwrap();
    write_pair(&directory.path().join("inputs/b"), ">seq\nATGAAATAA\n");
    write_pair(&directory.path().join("inputs/a"), ">seq\nATGGCTTAA\n");
    let output = analyse(directory.path(), &["inputs", "./inputs/a/sample.gff3"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let session = session(directory.path());
    assert_eq!(session.analysis.clusters.len(), 2);
    assert_eq!(
        session.analysis.clusters[0].loci[0].genes[0].translation,
        "MA"
    );
    assert_eq!(
        session.analysis.clusters[1].loci[0].genes[0].translation,
        "MK"
    );
}

#[cfg(unix)]
#[test]
fn discovery_skips_directory_symlink_cycles() {
    let directory = tempfile::tempdir().unwrap();
    write_pair(&directory.path().join("inputs"), ">seq\nATGGCTTAA\n");
    std::os::unix::fs::symlink(".", directory.path().join("inputs/loop")).unwrap();
    let output = analyse(directory.path(), &["inputs"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(session(directory.path()).analysis.clusters.len(), 1);
}

#[cfg(unix)]
#[test]
fn a_fasta_symlink_keeps_its_pairing_filename() {
    let directory = tempfile::tempdir().unwrap();
    fs::write(directory.path().join("sample.gff3"), GFF).unwrap();
    fs::write(directory.path().join("sequence.txt"), ">seq\nATGGCTTAA\n").unwrap();
    std::os::unix::fs::symlink("sequence.txt", directory.path().join("sample.fna")).unwrap();
    let output = analyse(directory.path(), &["sample.gff3", "sample.fna"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        session(directory.path()).analysis.clusters[0].loci[0].genes[0].translation,
        "MA"
    );
}
