use std::{fs, path::PathBuf};

use clap::Parser;
use clinker_core::parse_genbank;

/// Inspect GenBank input with the in-progress Rust implementation of clinker.
#[derive(Debug, Parser)]
#[command(version, about)]
struct Args {
    /// GenBank file to parse.
    input: PathBuf,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let bytes = fs::read(&args.input)?;
    let name = args.input.to_string_lossy();
    let cluster = parse_genbank(&name, &bytes)?;

    println!(
        "{}: {} locus/loci, {} CDS genes",
        cluster.name,
        cluster.loci.len(),
        cluster
            .loci
            .iter()
            .map(|locus| locus.genes.len())
            .sum::<usize>(),
    );

    Ok(())
}
