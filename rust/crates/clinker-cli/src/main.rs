use std::{fs, path::PathBuf};

use clap::Parser;
use clinker_core::{AnalysisOptions, DEFAULT_CONTIGUITY_WEIGHT, InputFile, analyse_genbank};

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

    /// Write the normalized synteny distance matrix as CSV.
    #[arg(long, value_name = "PATH")]
    matrix_out: Option<PathBuf>,

    /// Keep the input cluster order in plot JSON instead of synteny ordering.
    #[arg(long)]
    use_file_order: bool,
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
        let order = if args.use_file_order {
            (0..analysis.clusters.len()).collect()
        } else {
            analysis.cluster_order(DEFAULT_CONTIGUITY_WEIGHT)
        };
        let json = serde_json::to_string_pretty(&analysis.to_auto_arranged_plot_data(&order))?;
        fs::write(path, json)?;
    }

    if let Some(path) = args.matrix_out {
        fs::write(path, format_distance_matrix(&analysis))?;
    }

    Ok(())
}

fn format_distance_matrix(analysis: &clinker_core::Analysis) -> String {
    let matrix = analysis.synteny_distance_matrix(DEFAULT_CONTIGUITY_WEIGHT);
    let mut rows = vec![
        std::iter::once(String::new())
            .chain(analysis.clusters.iter().map(|cluster| cluster.name.clone()))
            .collect::<Vec<_>>(),
    ];
    rows.extend(matrix.iter().enumerate().map(|(index, row)| {
        std::iter::once(analysis.clusters[index].name.clone())
            .chain(row.iter().map(ToString::to_string))
            .collect()
    }));
    rows.into_iter()
        .map(|row| row.join(","))
        .collect::<Vec<_>>()
        .join("\n")
}
