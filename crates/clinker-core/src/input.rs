//! Dispatch named annotation bytes to the format-specific parsers.
//!
//! This module owns filename conventions and companion FASTA pairing. It does
//! not perform protein alignment; callers decide what to do with parsed clusters.

use std::collections::HashSet;

use crate::{
    AnalysisError, Cluster, InputFile, ParseOptions, parse_genbank_with_options, parse_gff3,
};

/// Parse inputs without comparing proteins.
///
/// Browser callers use this first to obtain cluster metadata, then distribute
/// independent pairwise comparisons to workers.
pub fn parse_input_files(files: &[InputFile<'_>]) -> Result<Vec<Cluster>, AnalysisError> {
    parse_input_files_with_options(files, &ParseOptions::default())
}

pub fn parse_input_files_with_options(
    files: &[InputFile<'_>],
    options: &ParseOptions,
) -> Result<Vec<Cluster>, AnalysisError> {
    parse_input_files_with_progress(files, options, |_| {})
}

/// Parse inputs with explicit options and report successfully processed files.
pub fn parse_input_files_with_progress<F>(
    files: &[InputFile<'_>],
    options: &ParseOptions,
    mut progress: F,
) -> Result<Vec<Cluster>, AnalysisError>
where
    F: FnMut(usize),
{
    let companion_indexes = files
        .iter()
        .filter(|file| is_gff(file.name))
        .filter_map(|gff| {
            let stem = file_stem(gff.name);
            files.iter().position(|candidate| {
                is_fasta(candidate.name) && file_stem(candidate.name).eq_ignore_ascii_case(stem)
            })
        })
        .collect::<HashSet<_>>();
    let mut clusters = Vec::with_capacity(files.len());
    for (index, file) in files.iter().enumerate() {
        if is_fasta(file.name) {
            if !companion_indexes.contains(&index) {
                return Err(AnalysisError::UnsupportedInput(format!(
                    "{} (no GFF3/GFF file with the same stem)",
                    file.name
                )));
            }
            progress(index + 1);
            continue;
        }
        if is_gff(file.name) {
            let stem = file_stem(file.name);
            let companion = files
                .iter()
                .find(|candidate| {
                    is_fasta(candidate.name) && file_stem(candidate.name).eq_ignore_ascii_case(stem)
                })
                .map(|candidate| candidate.bytes);
            clusters.push(
                parse_gff3(file.name, file.bytes, companion, options).map_err(|source| {
                    AnalysisError::ParseGff {
                        file_name: file.name.to_owned(),
                        source,
                    }
                })?,
            );
        } else if is_genbank(file.name) {
            clusters.extend(
                parse_genbank_with_options(file.name, file.bytes, options).map_err(|source| {
                    AnalysisError::Parse {
                        file_name: file.name.to_owned(),
                        source,
                    }
                })?,
            );
        } else {
            return Err(AnalysisError::UnsupportedInput(file.name.to_owned()));
        }
        progress(index + 1);
    }
    Ok(clusters)
}

fn extension(name: &str) -> Option<String> {
    name.rsplit_once('.')
        .map(|(_, extension)| extension.to_ascii_lowercase())
}

fn file_stem(name: &str) -> &str {
    name.rsplit_once('.').map_or(name, |(stem, _)| stem)
}

fn is_genbank(name: &str) -> bool {
    matches!(
        extension(name).as_deref(),
        Some("gb" | "gbk" | "genbank" | "gbf" | "gbff")
    )
}

fn is_gff(name: &str) -> bool {
    matches!(extension(name).as_deref(), Some("gff" | "gff3"))
}

fn is_fasta(name: &str) -> bool {
    matches!(
        extension(name).as_deref(),
        Some("fa" | "fsa" | "fna" | "fasta" | "faa")
    )
}
