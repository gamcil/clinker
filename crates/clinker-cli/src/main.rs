mod inputs;
mod outputs;
mod portable;
mod viewer;

use outputs::{Outputs, PlotFormat};
use portable::html as portable_html;

use std::{
    collections::BTreeMap,
    fs,
    io::{self, ErrorKind},
    path::{Path, PathBuf},
};

use clap::Parser;
use clinker_core::{
    AlignmentSettings, Analysis, AnalysisSession, ClinkerProject, GeneGroupDefinition,
    KmerPrefilter, ParseOptions, RowOrder, SequenceRange,
};

/// Analyse GenBank or GFF3 input with clinker.
#[derive(Debug, Parser)]
#[command(name = "clinker", version, about)]
struct Args {
    /// GenBank or GFF3 files/directories to analyse. Companion FASTA is discovered automatically.
    inputs: Vec<PathBuf>,

    /// Minimum protein identity required to retain a homology link (default: 0.30 for new sessions).
    #[arg(short, long, value_parser = parse_identity)]
    identity: Option<f32>,

    /// Number of alignment jobs to run in parallel (0 uses available CPUs).
    #[arg(short = 'j', long, default_value_t = 0)]
    jobs: usize,

    /// Skip protein pairs with too few shared k-mers; may miss distant homologues.
    #[arg(long, conflicts_with = "no_align")]
    prefilter: bool,

    /// Disable the saved approximate prefilter when reanalysing a session.
    #[arg(long, conflicts_with_all = ["prefilter", "no_align"])]
    no_prefilter: bool,

    /// Rebuild every comparison using saved translations and the requested settings.
    #[arg(long, requires = "session", conflicts_with = "open_project")]
    reanalyse: bool,

    /// Protein k-mer size for --prefilter (default: 3 for new sessions).
    #[arg(long, requires = "prefilter", value_parser = clap::value_parser!(u8).range(1..=8))]
    kmer_size: Option<u8>,

    /// Minimum distinct shared k-mers for --prefilter (default: 3 for new sessions).
    #[arg(long, requires = "prefilter", value_parser = clap::value_parser!(u32).range(1..))]
    min_shared_kmers: Option<u32>,

    /// Save analysis state here, or resume it when the path already exists.
    #[arg(short, long, value_name = "PATH")]
    session: Option<PathBuf>,

    /// Restrict a record to a 1-based inclusive range: RECORD:START-END.
    #[arg(short, long, value_name = "RECORD:START-END", num_args = 1..)]
    ranges: Vec<String>,

    /// Treat each GenBank record as a separate cluster.
    #[arg(long)]
    as_separate_clusters: bool,

    /// Keep the original origin of circular records.
    #[arg(long)]
    dont_set_origin: bool,

    /// Parse and plot genes without running protein alignments.
    #[arg(long)]
    no_align: bool,

    /// Two-column gene,function CSV used to seed plot groups.
    #[arg(long, value_name = "PATH")]
    gene_functions: Option<PathBuf>,

    /// Two-column function,colour CSV applied to function groups.
    #[arg(long, value_name = "PATH")]
    colour_map: Option<PathBuf>,

    /// Open the local viewer, or write portable HTML / plot JSON when PATH is given.
    #[arg(short, long, value_name = "PATH", num_args = 0..=1)]
    plot: Option<Option<PathBuf>>,

    /// Write a complete editable project for the web or local viewer.
    #[arg(long, value_name = "PATH")]
    project: Option<PathBuf>,

    /// Open an edited project in the local viewer without rerunning analysis.
    #[arg(long, value_name = "PATH", conflicts_with_all = ["inputs", "session", "plot", "project", "output", "matrix_out"])]
    open_project: Option<PathBuf>,

    /// Print the local viewer URL without launching a browser.
    #[arg(long)]
    no_browser: bool,

    /// Write the retained-link summary to a file instead of standard output.
    #[arg(short, long, value_name = "PATH")]
    output: Option<PathBuf>,

    /// Delimit summary columns instead of producing an aligned table.
    #[arg(short = 'l', long)]
    delimiter: Option<String>,

    /// Decimal places in identity and similarity output.
    #[arg(short = 'd', long, default_value_t = 2, value_parser = clap::value_parser!(u8).range(0..=12))]
    decimals: u8,

    /// Hide Query/Target/Identity/Similarity column headers.
    #[arg(long)]
    hide_link_headers: bool,

    /// Hide cluster-pair headings in the link summary.
    #[arg(long)]
    hide_alignment_headers: bool,

    /// Write the normalized synteny distance matrix as CSV.
    #[arg(long, value_name = "PATH")]
    matrix_out: Option<PathBuf>,

    /// Keep the input cluster order in plot JSON instead of synteny ordering.
    #[arg(long)]
    use_file_order: bool,

    /// Overwrite existing output files.
    #[arg(short, long)]
    force: bool,
}

fn alignment_settings(args: &Args, saved: Option<AlignmentSettings>) -> AlignmentSettings {
    if args.no_align {
        return AlignmentSettings::Skipped;
    }
    if !args.reanalyse && !args.prefilter && !args.no_prefilter {
        return saved.unwrap_or(AlignmentSettings::Aligned { prefilter: None });
    }
    let inherited = match saved {
        Some(AlignmentSettings::Aligned { prefilter }) => prefilter,
        _ => None,
    };
    let prefilter = if args.no_prefilter {
        None
    } else if args.prefilter {
        let defaults = inherited.unwrap_or_default();
        Some(KmerPrefilter {
            kmer_size: args
                .kmer_size
                .map(usize::from)
                .unwrap_or(defaults.kmer_size),
            min_shared_kmers: args
                .min_shared_kmers
                .map(|value| value as usize)
                .unwrap_or(defaults.min_shared_kmers),
        })
    } else {
        inherited
    };
    AlignmentSettings::Aligned { prefilter }
}

fn main() -> std::process::ExitCode {
    match run() {
        Ok(()) => std::process::ExitCode::SUCCESS,
        Err(error) => {
            eprintln!("Error: {error}");
            std::process::ExitCode::FAILURE
        }
    }
}

fn run() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    if let Some(path) = args.open_project.as_ref() {
        let project = ClinkerProject::from_json(&fs::read(path)?)?;
        let similarities = project.analysis.analysis.cluster_similarity_matrix();
        return viewer::serve(&project, &similarities, args.no_browser);
    }

    let session_exists = args.session.as_ref().is_some_and(|path| path.exists());
    if args.inputs.is_empty() && !session_exists {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            "provide at least one GenBank/GFF3 file, directory, or an existing --session",
        )
        .into());
    }
    if args.reanalyse && !session_exists {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            "--reanalyse requires an existing session",
        )
        .into());
    }
    let input_sources = inputs::discover(&args.inputs)?;
    if !args.inputs.is_empty() && input_sources.is_empty() {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "no supported GenBank or GFF3 files were found",
        )
        .into());
    }
    let mut outputs = Outputs::prepare(&args, &input_sources, session_exists)?;
    let parse_options = ParseOptions {
        ranges: parse_ranges(&args.ranges)?,
        as_separate_clusters: args.as_separate_clusters,
        set_origin: !args.dont_set_origin,
    };
    let clusters = inputs::parse(&input_sources, &parse_options)?;
    let supplied_group_definitions =
        read_group_definitions(args.gene_functions.as_deref(), args.colour_map.as_deref())?;
    let mut session = if session_exists {
        AnalysisSession::from_json(&fs::read(
            args.session.as_ref().expect("existing session path"),
        )?)?
    } else {
        AnalysisSession::new(
            Analysis {
                clusters: Vec::new(),
                links: Vec::new(),
            },
            args.identity.unwrap_or(0.30),
            alignment_settings(&args, None),
            Vec::new(),
        )?
    };
    let identity_cutoff = args.identity.unwrap_or(session.identity_cutoff());
    let alignment = alignment_settings(&args, session.alignment());
    alignment.validate(identity_cutoff)?;
    let explicit_settings =
        args.identity.is_some() || args.prefilter || args.no_prefilter || args.no_align;
    if session_exists && !args.reanalyse && (!clusters.is_empty() || explicit_settings) {
        session.ensure_settings(identity_cutoff, alignment)?;
    }
    if args.reanalyse {
        session.analysis.append_clusters_without_alignment(clusters);
        if let Some(options) = alignment.options(identity_cutoff) {
            session
                .analysis
                .reanalyse_with_options_and_jobs(options, args.jobs)?;
        } else {
            session.analysis.links.clear();
        }
        session = AnalysisSession::new(
            session.analysis,
            identity_cutoff,
            alignment,
            session.group_definitions,
        )?;
    } else if !clusters.is_empty() {
        if let Some(options) = alignment.options(identity_cutoff) {
            session
                .analysis
                .append_clusters_with_options_and_jobs(clusters, options, args.jobs)?;
        } else {
            session.analysis.append_clusters_without_alignment(clusters);
        }
    }
    if !session_exists || args.gene_functions.is_some() {
        session.group_definitions = supplied_group_definitions;
    }
    session.validate()?;
    let identity_cutoff = session.identity_cutoff();
    let analysis = &session.analysis;
    let group_definitions = &session.group_definitions;

    let summary = analysis.format_link_summary_with_options(
        usize::from(args.decimals),
        args.delimiter.as_deref(),
        !args.hide_alignment_headers,
        !args.hide_link_headers,
    );
    if let Some(output) = &mut outputs.summary {
        output.write(summary.as_bytes())?;
    }

    let figure = if args.plot.is_some() || args.project.is_some() {
        let row_order =
            if args.use_file_order || session.alignment() == Some(AlignmentSettings::Skipped) {
                RowOrder::Input
            } else {
                RowOrder::Similarity
            };
        Some(analysis.prepare_figure(row_order, group_definitions, |_| {})?)
    } else {
        None
    };
    let mut local_project = None;
    if let Some(figure) = &figure {
        let data = &figure.plot_data;
        if let Some(output) = &mut outputs.plot {
            let json = serde_json::to_string_pretty(data)?;
            let contents = match outputs
                .plot_format
                .as_ref()
                .expect("plot filename was validated")
            {
                PlotFormat::Json => json,
                PlotFormat::Html => portable_html(&json),
            };
            output.write(contents.as_bytes())?;
        }
        if outputs.project.is_some() || matches!(args.plot, Some(None)) {
            let project = ClinkerProject::new(session.clone(), data)?;
            if let Some(output) = &mut outputs.project {
                output.write(&project.to_json()?)?;
            }
            if matches!(args.plot, Some(None)) {
                local_project = Some((project, &figure.similarity_matrix));
            }
        }
    }

    if let Some(output) = &mut outputs.matrix {
        let matrix;
        let similarities = if let Some(figure) = &figure {
            &figure.similarity_matrix
        } else {
            matrix = analysis.cluster_similarity_matrix();
            &matrix
        };
        output.write(&outputs::distance_matrix(analysis, similarities)?)?;
    }
    outputs.commit(&session)?;
    if args.output.is_none() {
        if session.alignment() == Some(AlignmentSettings::Skipped) {
            println!(
                "Alignment disabled; parsed {} clusters",
                analysis.clusters.len()
            );
        } else if summary.is_empty() {
            println!("No retained links at identity cutoff {identity_cutoff:.2}");
        } else {
            println!("{summary}");
        }
    }

    if let Some((project, similarities)) = local_project {
        viewer::serve(&project, similarities, args.no_browser)?;
    }

    Ok(())
}

fn parse_identity(value: &str) -> Result<f32, String> {
    let identity = value
        .parse::<f32>()
        .map_err(|_| "identity must be a number between 0 and 1".to_owned())?;
    if identity.is_finite() && (0.0..=1.0).contains(&identity) {
        Ok(identity)
    } else {
        Err("identity must be a number between 0 and 1".to_owned())
    }
}

fn parse_ranges(values: &[String]) -> Result<BTreeMap<String, SequenceRange>, io::Error> {
    values
        .iter()
        .map(|value| {
            let (record, coordinates) = value.rsplit_once(':').ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("invalid range {value:?}; expected RECORD:START-END"),
                )
            })?;
            let (start, end) = coordinates.split_once('-').ok_or_else(|| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("invalid range {value:?}; expected RECORD:START-END"),
                )
            })?;
            let start = start.parse::<u64>().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("range start in {value:?} must be a positive integer"),
                )
            })?;
            let end = end.parse::<u64>().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("range end in {value:?} must be a positive integer"),
                )
            })?;
            if record.is_empty() || start == 0 || end < start {
                return Err(io::Error::new(
                    ErrorKind::InvalidInput,
                    format!(
                        "invalid range {value:?}; expected RECORD:START-END with 1 <= START <= END"
                    ),
                ));
            }
            Ok((
                record.to_owned(),
                SequenceRange {
                    start: start - 1,
                    end,
                },
            ))
        })
        .collect()
}

fn read_group_definitions(
    functions_path: Option<&Path>,
    colours_path: Option<&Path>,
) -> Result<Vec<GeneGroupDefinition>, Box<dyn std::error::Error>> {
    if functions_path.is_none() && colours_path.is_some() {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            "--colour-map requires --gene-functions",
        )
        .into());
    }
    let functions = functions_path
        .map(read_csv_pairs)
        .transpose()?
        .unwrap_or_default();
    let colours = colours_path
        .map(read_csv_pairs)
        .transpose()?
        .unwrap_or_default()
        .into_iter()
        .collect::<BTreeMap<_, _>>();
    Ok(group_definitions(functions, &colours))
}

fn read_csv_pairs(path: &Path) -> Result<Vec<(String, String)>, Box<dyn std::error::Error>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .trim(csv::Trim::All)
        .from_path(path)?;
    let mut pairs = Vec::new();
    for (index, row) in reader.records().enumerate() {
        let row = row?;
        if row.len() != 2 || row[0].is_empty() || row[1].is_empty() {
            return Err(io::Error::new(
                ErrorKind::InvalidData,
                format!(
                    "{} row {} must contain exactly two non-empty columns",
                    path.display(),
                    index + 1
                ),
            )
            .into());
        }
        pairs.push((row[0].to_owned(), row[1].to_owned()));
    }
    Ok(pairs)
}

fn group_definitions(
    functions: Vec<(String, String)>,
    colours: &BTreeMap<String, String>,
) -> Vec<GeneGroupDefinition> {
    let mut grouped = BTreeMap::<String, Vec<String>>::new();
    for (gene, function) in functions {
        grouped.entry(function).or_default().push(gene);
    }
    grouped
        .into_iter()
        .map(|(label, genes)| GeneGroupDefinition {
            colour: colours.get(&label).cloned(),
            label,
            genes,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{group_definitions, parse_identity, parse_ranges, portable_html};

    #[test]
    fn validates_prefilter_arguments() {
        use clap::Parser;
        assert!(
            super::Args::try_parse_from(["clinker", "--prefilter"])
                .unwrap()
                .prefilter
        );
        for args in [
            vec!["--kmer-size", "5"],
            vec!["--min-shared-kmers", "8"],
            vec!["--prefilter", "--kmer-size", "0"],
            vec!["--prefilter", "--kmer-size", "9"],
            vec!["--prefilter", "--min-shared-kmers", "0"],
            vec!["--prefilter", "--no-align"],
        ] {
            assert!(super::Args::try_parse_from(std::iter::once("clinker").chain(args)).is_err());
        }
    }

    #[test]
    fn validates_identity_range() {
        assert_eq!(parse_identity("0.3"), Ok(0.3));
        assert!(parse_identity("-0.1").is_err());
        assert!(parse_identity("1.1").is_err());
        assert!(parse_identity("NaN").is_err());
    }

    #[test]
    fn portable_html_embeds_data_without_a_script_terminator() {
        let html = portable_html(r#"{"label":"</script>"}"#);

        assert!(html.starts_with("<!doctype html>"));
        assert!(html.contains("ClusterMap.ClusterMap()"));
        assert!(html.contains("#plot{position:relative;width:100%;height:100%;overflow:hidden}"));
        assert!(!html.contains("min-height:100vh"));
        assert!(html.contains(r#"\u003c/script\u003e"#));
        assert!(!html.contains(r#"{"label":"</script>"}"#));
    }

    #[test]
    fn converts_cli_ranges_to_internal_coordinates() {
        let ranges = parse_ranges(&["record:5-12".to_owned()]).unwrap();
        let range = ranges["record"];

        assert_eq!((range.start, range.end), (4, 12));
        assert!(parse_ranges(&["record:0-12".to_owned()]).is_err());
        assert!(parse_ranges(&["record:12-5".to_owned()]).is_err());
    }

    #[test]
    fn combines_function_rows_and_colours_into_group_definitions() {
        let definitions = group_definitions(
            vec![
                ("gene-a".into(), "Enzyme".into()),
                ("gene-b".into(), "Enzyme".into()),
            ],
            &BTreeMap::from([("Enzyme".into(), "#123456".into())]),
        );

        assert_eq!(definitions.len(), 1);
        assert_eq!(definitions[0].label, "Enzyme");
        assert_eq!(definitions[0].genes, ["gene-a", "gene-b"]);
        assert_eq!(definitions[0].colour.as_deref(), Some("#123456"));
    }
}
