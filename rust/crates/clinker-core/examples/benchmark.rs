use std::{
    error::Error,
    fs,
    io::{self, Read},
    time::Instant,
};

use clinker_core::{Analysis, AnalysisOptions, InputFile, KmerPrefilter, parse_input_files};
use serde::Deserialize;
use serde_json::json;

#[derive(Deserialize)]
struct Request {
    files: Vec<String>,
    workers: usize,
    identity: f32,
    prefilter: bool,
    kmer_size: usize,
    min_shared_kmers: usize,
    inspect: bool,
}

fn main() -> Result<(), Box<dyn Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let request: Request = serde_json::from_str(&input)?;
    if request.workers == 0 {
        return Err("workers must be positive".into());
    }
    let staged = request
        .files
        .iter()
        .map(|path| {
            Ok((
                std::path::Path::new(path)
                    .file_name()
                    .ok_or("invalid filename")?
                    .to_string_lossy()
                    .into_owned(),
                fs::read(path)?,
            ))
        })
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let inputs = staged
        .iter()
        .map(|(name, bytes)| InputFile { name, bytes })
        .collect::<Vec<_>>();
    let start = Instant::now();
    let clusters = parse_input_files(&inputs)?;
    let parse_ms = start.elapsed().as_secs_f64() * 1000.0;
    if request.inspect {
        let proteins = clusters
            .iter()
            .map(|cluster| {
                cluster
                    .loci
                    .iter()
                    .flat_map(|locus| locus.genes.iter().map(|gene| gene.translation.clone()))
                    .collect::<Vec<_>>()
            })
            .collect::<Vec<_>>();
        println!("{}", json!({"proteins": proteins}));
        return Ok(());
    }
    let options = AnalysisOptions {
        identity_cutoff: request.identity,
        prefilter: request.prefilter.then_some(KmerPrefilter {
            kmer_size: request.kmer_size,
            min_shared_kmers: request.min_shared_kmers,
        }),
    };
    let mut analysis = Analysis {
        clusters: Vec::new(),
        links: Vec::new(),
    };
    let stats =
        analysis.append_clusters_with_options_and_jobs(clusters, options, request.workers)?;
    let total_ms = start.elapsed().as_secs_f64() * 1000.0;
    // Work counts are diagnostics, computed after the elapsed-time boundary.
    let clusters = &analysis.clusters;
    let proteins = clusters
        .iter()
        .flat_map(|cluster| {
            cluster
                .loci
                .iter()
                .flat_map(|locus| locus.genes.iter().map(|gene| gene.translation.as_bytes()))
        })
        .collect::<Vec<_>>();
    let protein_count = proteins.len();
    let amino_acids = proteins.iter().map(|p| p.len()).sum::<usize>();
    println!(
        "{}",
        json!({"parse_ms":parse_ms,"prepare_ms":null,
        "prefilter_ms":null,"alignment_ms":null,"total_ms":total_ms,
        "proteins":protein_count,"amino_acids":amino_acids,"clusters":clusters.len(),
        "possible_pairs":stats.counts.possible_pairs,"eligible_pairs":null,"candidate_pairs":stats.counts.planned_pairs,
        "aligned_pairs":stats.aligned_pairs,"retained_links":analysis.links.len()})
    );
    Ok(())
}
