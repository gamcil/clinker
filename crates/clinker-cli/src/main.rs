use std::{
    collections::BTreeMap,
    fs,
    io::{self, ErrorKind},
    path::{Path, PathBuf},
};

use clap::Parser;
use clinker_core::{
    Analysis, AnalysisOptions, AnalysisSession, ClinkerProject, GeneGroupDefinition, InputFile,
    KmerPrefilter, ParseOptions, SequenceRange, analyse_genbank_with_options_and_jobs,
    parse_input_files_with_options,
};

/// Analyse GenBank or GFF3 input with clinker.
#[derive(Debug, Parser)]
#[command(version, about)]
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

    /// Protein k-mer size for --prefilter.
    #[arg(long, default_value_t = 3, requires = "prefilter", value_parser = clap::value_parser!(u8).range(1..=8))]
    kmer_size: u8,

    /// Minimum distinct shared k-mers for --prefilter.
    #[arg(long, default_value_t = 3, requires = "prefilter", value_parser = clap::value_parser!(u32).range(1..))]
    min_shared_kmers: u32,

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

    /// Write portable HTML (.html) or clustermap.js data (.json).
    #[arg(short, long, value_name = "PATH")]
    plot: Option<PathBuf>,

    /// Write a complete editable project for the web or local viewer.
    #[arg(long, value_name = "PATH")]
    project: Option<PathBuf>,

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

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();
    let session_exists = args.session.as_ref().is_some_and(|path| path.exists());
    if args.inputs.is_empty() && !session_exists {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            "provide at least one GenBank/GFF3 file, directory, or an existing --session",
        )
        .into());
    }
    let input_paths = discover_input_paths(&args.inputs)?;
    if !args.inputs.is_empty() && input_paths.is_empty() {
        return Err(io::Error::new(
            ErrorKind::NotFound,
            "no supported GenBank or GFF3 files were found",
        )
        .into());
    }
    let owned_inputs = input_paths
        .iter()
        .map(|path| {
            let name = path
                .file_name()
                .unwrap_or(path.as_os_str())
                .to_string_lossy()
                .into_owned();
            Ok((name, fs::read(path)?))
        })
        .collect::<Result<Vec<_>, std::io::Error>>()?;
    let inputs = owned_inputs
        .iter()
        .map(|(name, bytes)| InputFile { name, bytes })
        .collect::<Vec<_>>();
    let parse_options = ParseOptions {
        ranges: parse_ranges(&args.ranges)?,
        as_separate_clusters: args.as_separate_clusters,
        set_origin: !args.dont_set_origin,
    };
    let supplied_group_definitions =
        read_group_definitions(args.gene_functions.as_deref(), args.colour_map.as_deref())?;
    let replace_group_definitions = args.gene_functions.is_some();
    let prefilter = args.prefilter.then_some(KmerPrefilter {
        kmer_size: usize::from(args.kmer_size),
        min_shared_kmers: args.min_shared_kmers as usize,
    });
    let mut session_was_modified = false;
    let session = if session_exists {
        let path = args
            .session
            .as_ref()
            .expect("checked existing session path");
        let mut session = AnalysisSession::from_json(&fs::read(path)?)?;
        if !inputs.is_empty() {
            let identity_cutoff = args.identity.unwrap_or(session.identity_cutoff);
            let clusters = parse_input_files_with_options(&inputs, &parse_options)?;
            if args.no_align {
                session.analysis.append_clusters_without_alignment(clusters);
            } else {
                session.analysis.append_clusters_with_options_and_jobs(
                    clusters,
                    AnalysisOptions {
                        identity_cutoff,
                        prefilter,
                    },
                    args.jobs,
                )?;
            }
            session.identity_cutoff = identity_cutoff;
            session_was_modified = true;
        }
        if replace_group_definitions {
            session.group_definitions = supplied_group_definitions;
            session_was_modified = true;
        }
        session
    } else {
        let identity_cutoff = args.identity.unwrap_or(0.30);
        let analysis = if args.no_align {
            Analysis {
                clusters: parse_input_files_with_options(&inputs, &parse_options)?,
                links: Vec::new(),
            }
        } else {
            analyse_genbank_with_options_and_jobs(
                &inputs,
                AnalysisOptions {
                    identity_cutoff,
                    prefilter,
                },
                &parse_options,
                args.jobs,
            )?
        };
        session_was_modified = args.session.is_some();
        AnalysisSession::new(analysis, identity_cutoff, supplied_group_definitions)?
    };
    session.validate()?;
    if session_was_modified && let Some(path) = args.session.as_ref() {
        fs::write(path, session.to_json()?)?;
    }
    let identity_cutoff = session.identity_cutoff;
    let analysis = &session.analysis;
    let group_definitions = &session.group_definitions;

    let summary = analysis.format_link_summary_with_options(
        usize::from(args.decimals),
        args.delimiter.as_deref(),
        !args.hide_alignment_headers,
        !args.hide_link_headers,
    );
    if let Some(path) = args.output.as_ref() {
        write_output(path, summary.as_bytes(), args.force)?;
    } else if args.no_align {
        println!(
            "Alignment disabled; parsed {} clusters",
            analysis.clusters.len()
        );
    } else if summary.is_empty() {
        println!("No retained links at identity cutoff {identity_cutoff:.2}");
    } else {
        println!("{summary}");
    }

    if args.plot.is_some() || args.project.is_some() {
        let order = if args.use_file_order || args.no_align {
            (0..analysis.clusters.len()).collect()
        } else {
            analysis.cluster_order()
        };
        let data =
            analysis.to_auto_arranged_plot_data_with_group_definitions(&order, group_definitions);
        if let Some(path) = args.plot.as_ref() {
            let json = serde_json::to_string_pretty(&data)?;
            let output = match path
                .extension()
                .and_then(|extension| extension.to_str())
                .map(str::to_ascii_lowercase)
                .as_deref()
            {
                Some("json") => json,
                Some("html") | Some("htm") => portable_html(&json),
                _ => {
                    return Err(io::Error::new(
                        ErrorKind::InvalidInput,
                        "--plot output must end in .json, .html, or .htm",
                    )
                    .into());
                }
            };
            write_output(path, output.as_bytes(), args.force)?;
        }
        if let Some(path) = args.project.as_ref() {
            let project = ClinkerProject::new(session.clone(), &data)?;
            write_output(path, &project.to_json()?, args.force)?;
        }
    }

    if let Some(path) = args.matrix_out.as_ref() {
        write_output(
            path,
            format_distance_matrix(analysis).as_bytes(),
            args.force,
        )?;
    }

    Ok(())
}

fn discover_input_paths(inputs: &[PathBuf]) -> io::Result<Vec<PathBuf>> {
    let mut primary = Vec::new();
    for input in inputs {
        collect_input_path(input, true, &mut primary)?;
    }
    let mut discovered = Vec::new();
    for path in primary {
        push_unique(&mut discovered, path.clone());
        if is_gff_path(&path)
            && let Some(companion) = fasta_companion(&path)?
        {
            push_unique(&mut discovered, companion);
        }
    }
    Ok(discovered)
}

fn collect_input_path(path: &Path, explicit: bool, output: &mut Vec<PathBuf>) -> io::Result<()> {
    if path.is_dir() {
        let mut entries = fs::read_dir(path)?
            .map(|entry| entry.map(|entry| entry.path()))
            .collect::<Result<Vec<_>, _>>()?;
        entries.sort();
        for entry in entries {
            collect_input_path(&entry, false, output)?;
        }
    } else if is_primary_input(path) || (explicit && is_fasta_path(path)) {
        output.push(path.to_owned());
    } else if explicit {
        return Err(io::Error::new(
            ErrorKind::InvalidInput,
            format!(
                "unsupported input {}; expected GenBank, GFF3/GFF, or FASTA",
                path.display()
            ),
        ));
    }
    Ok(())
}

fn extension(path: &Path) -> Option<String> {
    path.extension()
        .and_then(|extension| extension.to_str())
        .map(str::to_ascii_lowercase)
}

fn is_primary_input(path: &Path) -> bool {
    matches!(
        extension(path).as_deref(),
        Some("gb" | "gbk" | "genbank" | "gbf" | "gbff" | "gff" | "gff3")
    )
}

fn is_gff_path(path: &Path) -> bool {
    matches!(extension(path).as_deref(), Some("gff" | "gff3"))
}

fn is_fasta_path(path: &Path) -> bool {
    matches!(
        extension(path).as_deref(),
        Some("fa" | "fsa" | "fna" | "fasta" | "faa")
    )
}

fn fasta_companion(gff: &Path) -> io::Result<Option<PathBuf>> {
    let parent = gff.parent().unwrap_or_else(|| Path::new("."));
    let stem = gff.file_stem().and_then(|stem| stem.to_str()).unwrap_or("");
    let mut entries = fs::read_dir(parent)?
        .map(|entry| entry.map(|entry| entry.path()))
        .collect::<Result<Vec<_>, _>>()?;
    entries.sort();
    Ok(entries.into_iter().find(|candidate| {
        is_fasta_path(candidate)
            && candidate
                .file_stem()
                .and_then(|candidate_stem| candidate_stem.to_str())
                .is_some_and(|candidate_stem| candidate_stem.eq_ignore_ascii_case(stem))
    }))
}

fn push_unique(paths: &mut Vec<PathBuf>, path: PathBuf) {
    if !paths.contains(&path) {
        paths.push(path);
    }
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
            let start = start.parse::<usize>().map_err(|_| {
                io::Error::new(
                    ErrorKind::InvalidInput,
                    format!("range start in {value:?} must be a positive integer"),
                )
            })?;
            let end = end.parse::<usize>().map_err(|_| {
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

fn write_output(path: &Path, contents: &[u8], force: bool) -> io::Result<()> {
    if !force && path.exists() {
        return Err(io::Error::new(
            ErrorKind::AlreadyExists,
            format!(
                "{} already exists; pass --force to overwrite it",
                path.display()
            ),
        ));
    }
    fs::write(path, contents)
}

fn portable_html(plot_json: &str) -> String {
    const D3: &str = include_str!("../assets/d3.min.js");
    const CLUSTERMAP: &str = include_str!("../assets/clustermap.min.js");
    let safe_json = plot_json
        .replace('&', "\\u0026")
        .replace('<', "\\u003c")
        .replace('>', "\\u003e")
        .replace('\u{2028}', "\\u2028")
        .replace('\u{2029}', "\\u2029");
    let mut html = String::with_capacity(D3.len() + CLUSTERMAP.len() + safe_json.len() + 1_024);
    html.push_str(
        "<!doctype html><html lang=\"en\"><head><meta charset=\"utf-8\">\
         <meta name=\"viewport\" content=\"width=device-width,initial-scale=1\">\
         <title>clinker plot</title><style>html,body{width:100%;height:100%;margin:0;overflow:hidden;font-family:system-ui,sans-serif}\
         #plot{position:relative;width:100%;height:100%;overflow:hidden}svg.clusterMap{display:block}</style></head><body>\
         <div id=\"plot\"></div><script>",
    );
    html.push_str(D3);
    html.push_str("</script><script>");
    html.push_str(CLUSTERMAP);
    html.push_str("</script><script>const data=");
    html.push_str(&safe_json);
    html.push_str(
        ";const chart=ClusterMap.ClusterMap();d3.select('#plot').datum(data).call(chart);</script></body></html>",
    );
    html
}

fn format_distance_matrix(analysis: &clinker_core::Analysis) -> String {
    let matrix = analysis.synteny_distance_matrix();
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

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::{
        discover_input_paths, group_definitions, parse_identity, parse_ranges, portable_html,
    };

    #[test]
    fn validates_prefilter_arguments() {
        use clap::Parser;
        assert!(
            super::Args::try_parse_from(["clinker-rs", "--prefilter"])
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
            assert!(
                super::Args::try_parse_from(std::iter::once("clinker-rs").chain(args)).is_err()
            );
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

    #[test]
    fn discovers_gff_companion_fasta() {
        let fixture = std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../clinker-core/tests/fixtures/features.gff3");
        let paths = discover_input_paths(&[fixture]).unwrap();

        assert_eq!(paths.len(), 2);
        assert_eq!(
            paths[0].extension().and_then(|value| value.to_str()),
            Some("gff3")
        );
        assert_eq!(
            paths[1].extension().and_then(|value| value.to_str()),
            Some("fna")
        );

        let directory = paths[0].parent().unwrap().to_owned();
        let discovered = discover_input_paths(&[directory]).unwrap();
        for path in &paths {
            assert_eq!(discovered.iter().filter(|found| *found == path).count(), 1);
        }
    }
}
