use std::{fs, process::Command};

#[test]
fn discovers_the_companion_fasta_and_writes_gff_plot_data() {
    let manifest = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let gff = manifest.join("tests/fixtures/features.gff3");
    let output = std::env::temp_dir().join(format!(
        "clinker-gff-{}-{}.json",
        std::process::id(),
        std::thread::current().name().unwrap_or("test")
    ));
    let result = Command::new(env!("CARGO_BIN_EXE_clinker-rs"))
        .arg(gff)
        .arg("--no-align")
        .arg("--plot")
        .arg(&output)
        .output()
        .expect("run clinker CLI");

    assert!(
        result.status.success(),
        "stderr: {}",
        String::from_utf8_lossy(&result.stderr)
    );
    let plot: serde_json::Value =
        serde_json::from_slice(&fs::read(&output).expect("plot data was written")).unwrap();
    assert_eq!(plot["clusters"][0]["loci"].as_array().unwrap().len(), 2);
    fs::remove_file(output).unwrap();
}
