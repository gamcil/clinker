use std::{env, error::Error, fs, path::Path, time::Instant};

use clinker_core::{
    AnalysisOptions, InputFile, KmerPrefilter, ParseOptions, RowOrder, analyse_with_jobs,
};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<String> = env::args().collect();
    if args.len() != 9 {
        return Err("expected fixture, mode, identity, threads, kmer size, shared kmers, figure path, metrics path".into());
    }
    let mut paths = fs::read_dir(&args[1])?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<Result<Vec<_>, _>>()?;
    paths.retain(|path| path.extension().is_some_and(|ext| ext == "gbk"));
    paths.sort();
    if paths.is_empty() {
        return Err("no GenBank files found".into());
    }
    let prefilter = if args[2] == "rust-kmer" {
        Some(KmerPrefilter {
            kmer_size: args[5].parse()?,
            min_shared_kmers: args[6].parse()?,
        })
    } else {
        None
    };
    let options = AnalysisOptions {
        identity_cutoff: args[3].parse()?,
        prefilter,
    };
    let threads = args[4].parse()?;

    let started = Instant::now();
    let contents = paths.iter().map(fs::read).collect::<Result<Vec<_>, _>>()?;
    let files = paths
        .iter()
        .zip(&contents)
        .map(|(path, bytes)| InputFile {
            name: path.file_name().unwrap().to_str().unwrap(),
            bytes,
        })
        .collect::<Vec<_>>();
    let analysis = analyse_with_jobs(&files, options, &ParseOptions::default(), threads)?;
    let figure = analysis.prepare_figure(RowOrder::Similarity, &[], |_| {})?;
    let json = serde_json::to_vec(&figure.plot_data)?;
    let elapsed = started.elapsed().as_secs_f64();

    fs::write(Path::new(&args[7]), &json)?;
    fs::write(
        Path::new(&args[8]),
        serde_json::to_vec(&serde_json::json!({
            "elapsed_seconds": elapsed, "retained_links": figure.plot_data.links.len(),
            "groups": figure.plot_data.groups.len(), "json_bytes": json.len()
        }))?,
    )?;
    Ok(())
}
