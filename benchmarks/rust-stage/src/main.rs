use std::{env, fs, hint::black_box, time::Instant};

use clinker_core::{
    Analysis, AnalysisOptions, InputFile, KmerPrefilter, ParseOptions,
    parse_input_files_with_options,
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut args = env::args().skip(1);
    let dataset = args.next().ok_or("missing fixture directory")?;
    let identity: f32 = args.next().ok_or("missing identity cutoff")?.parse()?;
    let jobs: usize = args.next().ok_or("missing thread count")?.parse()?;
    let filter = args.next().as_deref() == Some("kmer");
    let kmer_size: usize = args.next().ok_or("missing k-mer size")?.parse()?;
    let min_shared_kmers: usize = args.next().ok_or("missing shared k-mer count")?.parse()?;
    let output = args.next().ok_or("missing result path")?;

    let mut paths = fs::read_dir(dataset)?
        .filter_map(Result::ok)
        .map(|entry| entry.path())
        .filter(|path| path.extension().is_some_and(|extension| extension == "gbk"))
        .collect::<Vec<_>>();
    paths.sort();
    let bytes = paths
        .iter()
        .map(fs::read)
        .collect::<Result<Vec<_>, _>>()?;
    let names = paths
        .iter()
        .map(|path| path.to_str().ok_or("non-UTF8 input path"))
        .collect::<Result<Vec<_>, _>>()?;
    let inputs = paths
        .iter()
        .zip(names)
        .zip(&bytes)
        .map(|((_, name), bytes)| InputFile {
            name,
            bytes,
        })
        .collect::<Vec<_>>();
    let clusters = parse_input_files_with_options(&inputs, &ParseOptions::default())?;
    let options = AnalysisOptions {
        identity_cutoff: identity,
        prefilter: filter.then_some(KmerPrefilter {
            kmer_size,
            min_shared_kmers,
        }),
    };

    let mut analysis = Analysis {
        clusters: Vec::new(),
        links: Vec::new(),
    };
    let started = Instant::now();
    analysis.append_clusters_with_jobs(black_box(clusters), options, jobs)?;
    let elapsed = started.elapsed().as_secs_f64();
    fs::write(
        output,
        format!(
            "{}\n",
            serde_json::json!({
                "elapsed_seconds": elapsed,
                "retained_links": analysis.links.len(),
            })
        ),
    )?;
    Ok(())
}
