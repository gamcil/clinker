use std::{fs, path::PathBuf};

use clap::Parser;
use clinker_core::{AnalysisOptions, InputFile, analyse_genbank};

/// Inspect GenBank input with the in-progress Rust implementation of clinker.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Args {
    /// GenBank files to analyse.
    #[arg(required = true)]
    inputs: Vec<PathBuf>,

    /// Minimum protein identity required to retain a homology link.
    #[arg(short, long, default_value_t = 0.30)]
    identity: f32,

    /// Write data compatible with the existing clustermap.js renderer.
    #[arg(long, value_name = "PATH")]
    plot: Option<PathBuf>,

    /// Write the retained-link summary to a file instead of standard output.
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let owned_inputs = args
        .inputs
        .iter()
        .map(|path| Ok((path.to_string_lossy().into_owned(), fs::read(path)?)))
        .collect::<Result<Vec<_>, std::io::Error>>()?;
    let inputs = owned_inputs
        .iter()
        .map(|(name, bytes)| InputFile { name, bytes })
        .collect::<Vec<_>>();
    let analysis = analyse_genbank(
        &inputs,
        AnalysisOptions {
            identity_cutoff: args.identity,
        },
    )?;

    let summary = analysis.format_link_summary();
    if let Some(path) = args.output {
        fs::write(path, summary)?;
    } else if summary.is_empty() {
        println!("No retained links at identity cutoff {:.2}", args.identity);
    } else {
        println!("{summary}");
    }

    if let Some(path) = args.plot {
        let json = serde_json::to_string_pretty(&analysis.to_plot_data())?;
        fs::write(path, json)?;
    }

    Ok(())
}
