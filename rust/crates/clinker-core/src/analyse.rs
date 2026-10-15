use std::collections::HashSet;

use serde::{Deserialize, Serialize};
use thiserror::Error;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

use crate::align::ProteinAligner;
use crate::{Cluster, ParseOptions, ProteinMatch, parse_genbank_with_options, parse_gff3};

#[derive(Debug, Clone, Copy)]
pub struct InputFile<'a> {
    pub name: &'a str,
    pub bytes: &'a [u8],
}

#[derive(Debug, Clone, Copy)]
pub struct AnalysisOptions {
    pub identity_cutoff: f32,
}

impl Default for AnalysisOptions {
    fn default() -> Self {
        Self {
            identity_cutoff: 0.30,
        }
    }
}

impl Analysis {
    /// Format alignment output as in original clinker.
    /// For structured results use [`Self::links`] or [`Self::to_plot_data`].
    pub fn format_link_summary(&self) -> String {
        self.format_link_summary_with_options(4, None, true, true)
    }

    /// Format retained links with the same presentation controls as the
    /// Python command line interface.
    pub fn format_link_summary_with_options(
        &self,
        decimals: usize,
        delimiter: Option<&str>,
        alignment_headers: bool,
        link_headers: bool,
    ) -> String {
        let mut sections = Vec::new();
        for query_cluster in 0..self.clusters.len() {
            for target_cluster in query_cluster + 1..self.clusters.len() {
                let links = self
                    .links
                    .iter()
                    .filter(|link| {
                        link.query.cluster == query_cluster && link.target.cluster == target_cluster
                    })
                    .collect::<Vec<_>>();
                if links.is_empty() {
                    continue;
                }
                let header = format!(
                    "{} vs {}",
                    self.clusters[query_cluster].name, self.clusters[target_cluster].name
                );
                let mut rows = Vec::new();
                if link_headers {
                    rows.push([
                        "Query".to_owned(),
                        "Target".to_owned(),
                        "Identity".to_owned(),
                        "Similarity".to_owned(),
                    ]);
                }
                rows.extend(links.into_iter().map(|link| {
                    [
                        self.gene_label(link.query).to_owned(),
                        self.gene_label(link.target).to_owned(),
                        format!("{:.decimals$}", link.identity),
                        format!("{:.decimals$}", link.similarity),
                    ]
                }));
                let table = if let Some(delimiter) = delimiter {
                    rows.into_iter()
                        .map(|row| row.join(delimiter))
                        .collect::<Vec<_>>()
                        .join("\n")
                } else {
                    let widths = std::array::from_fn::<_, 4, _>(|column| {
                        rows.iter().map(|row| row[column].len()).max().unwrap_or(0)
                    });
                    rows.into_iter()
                        .map(|row| {
                            row.into_iter()
                                .enumerate()
                                .map(|(column, value)| {
                                    format!("{value:<width$}", width = widths[column])
                                })
                                .collect::<Vec<_>>()
                                .join("  ")
                        })
                        .collect::<Vec<_>>()
                        .join("\n")
                };

                sections.push(if alignment_headers {
                    format!("{header}\n{}\n{table}", "-".repeat(header.len()))
                } else {
                    table
                });
            }
        }

        sections.join("\n\n")
    }

    fn gene_label(&self, reference: GeneRef) -> &str {
        &self.clusters[reference.cluster].loci[reference.locus].genes[reference.gene].label
    }
}

#[derive(Debug, Error)]
pub enum AnalysisError {
    #[error("failed to parse {file_name}: {source}")]
    Parse {
        file_name: String,
        #[source]
        source: crate::ParseError,
    },
    #[error("failed to parse {file_name}: {source}")]
    ParseGff {
        file_name: String,
        #[source]
        source: crate::GffError,
    },
    #[error("unsupported input file {0}; expected GenBank, GFF3/GFF, or a companion FASTA")]
    UnsupportedInput(String),
    #[cfg(not(target_arch = "wasm32"))]
    #[error("failed to create the alignment thread pool: {0}")]
    ThreadPool(#[from] rayon::ThreadPoolBuildError),
}

/// A stable position within this analysis result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Deserialize, Serialize)]
pub struct GeneRef {
    pub cluster: usize,
    pub locus: usize,
    pub gene: usize,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Link {
    pub query: GeneRef,
    pub target: GeneRef,
    pub identity: f32,
    pub similarity: f32,
}

/// A retained match between two positions in a protein-comparison tile.
///
/// The indices refer to the caller-provided query and target slices rather
/// than loci or clusters. Browser workers use this compact result while the
/// coordinator retains the corresponding gene references.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProteinTileMatch {
    pub query_index: usize,
    pub target_index: usize,
    pub identity: f32,
    pub similarity: f32,
}

/// A pair of positions in caller-provided protein slices.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProteinPair {
    pub query_index: usize,
    pub target_index: usize,
}

#[derive(Debug, Clone, PartialEq, Deserialize, Serialize)]
pub struct Analysis {
    pub clusters: Vec<Cluster>,
    pub links: Vec<Link>,
}

impl Analysis {
    /// Append clusters and compare only pairs involving at least one newly
    /// appended cluster. Existing links and cluster indices remain stable.
    pub fn append_clusters(&mut self, clusters: Vec<Cluster>, identity_cutoff: f32) {
        let first_new_cluster = self.clusters.len();
        self.clusters.extend(clusters);
        self.links.extend(cross_cluster_links_from(
            &self.clusters,
            first_new_cluster,
            identity_cutoff,
        ));
    }

    /// Append clusters and compare new cluster pairs using a local native
    /// thread pool. `jobs == 0` uses the number of available CPU threads.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn append_clusters_with_jobs(
        &mut self,
        clusters: Vec<Cluster>,
        identity_cutoff: f32,
        jobs: usize,
    ) -> Result<(), AnalysisError> {
        let pool = alignment_pool(jobs)?;
        let first_new_cluster = self.clusters.len();
        self.clusters.extend(clusters);
        self.links.extend(cross_cluster_links_parallel_from(
            &self.clusters,
            first_new_cluster,
            identity_cutoff,
            &pool,
        ));
        Ok(())
    }

    /// Append parsed clusters without producing protein links.
    pub fn append_clusters_without_alignment(&mut self, clusters: Vec<Cluster>) {
        self.clusters.extend(clusters);
    }
}

/// Parse supported annotation inputs and retain cross-cluster gene links at the identity
/// cutoff.
pub fn analyse_genbank(
    files: &[InputFile<'_>],
    options: AnalysisOptions,
) -> Result<Analysis, AnalysisError> {
    analyse_genbank_with_options(files, options, &ParseOptions::default())
}

/// Analyse GenBank inputs with explicit parsing behavior.
pub fn analyse_genbank_with_options(
    files: &[InputFile<'_>],
    options: AnalysisOptions,
    parse_options: &ParseOptions,
) -> Result<Analysis, AnalysisError> {
    let clusters = parse_input_files_with_options(files, parse_options)?;

    Ok(Analysis {
        links: cross_cluster_links(&clusters, options.identity_cutoff),
        clusters,
    })
}

/// Analyse annotation inputs using a local native alignment thread pool.
///
/// `jobs == 0` uses the number of available CPU threads. Results retain the
/// same deterministic cluster-pair ordering as the serial implementation.
#[cfg(not(target_arch = "wasm32"))]
pub fn analyse_genbank_with_options_and_jobs(
    files: &[InputFile<'_>],
    options: AnalysisOptions,
    parse_options: &ParseOptions,
    jobs: usize,
) -> Result<Analysis, AnalysisError> {
    let clusters = parse_input_files_with_options(files, parse_options)?;
    let pool = alignment_pool(jobs)?;

    Ok(Analysis {
        links: cross_cluster_links_parallel_from(&clusters, 0, options.identity_cutoff, &pool),
        clusters,
    })
}

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
    parse_input_files_with_options_and_progress(files, options, |_| {})
}

/// Parse annotation inputs and report the number of successfully processed input files.
///
/// Browser workers use the callback to surface progress without moving file
/// parsing out of the shared core implementation.
pub fn parse_input_files_with_progress<F>(
    files: &[InputFile<'_>],
    progress: F,
) -> Result<Vec<Cluster>, AnalysisError>
where
    F: FnMut(usize),
{
    parse_input_files_with_options_and_progress(files, &ParseOptions::default(), progress)
}

pub fn parse_input_files_with_options_and_progress<F>(
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

fn cross_cluster_links(clusters: &[Cluster], identity_cutoff: f32) -> Vec<Link> {
    cross_cluster_links_from(clusters, 0, identity_cutoff)
}

fn cross_cluster_links_from(
    clusters: &[Cluster],
    first_new_cluster: usize,
    identity_cutoff: f32,
) -> Vec<Link> {
    let mut links = Vec::new();

    for query_cluster_index in 0..clusters.len() {
        let first_target = (query_cluster_index + 1).max(first_new_cluster);
        for target_cluster_index in first_target..clusters.len() {
            links.extend(
                analyse_cluster_pair(
                    &clusters[query_cluster_index],
                    &clusters[target_cluster_index],
                    identity_cutoff,
                )
                .into_iter()
                .map(|mut link| {
                    link.query.cluster = query_cluster_index;
                    link.target.cluster = target_cluster_index;
                    link
                }),
            );
        }
    }

    links
}

#[cfg(not(target_arch = "wasm32"))]
fn cross_cluster_links_parallel_from(
    clusters: &[Cluster],
    first_new_cluster: usize,
    identity_cutoff: f32,
    pool: &rayon::ThreadPool,
) -> Vec<Link> {
    let tasks = alignment_tasks(clusters, first_new_cluster);
    let max_protein_length = clusters
        .iter()
        .map(max_protein_length)
        .max()
        .unwrap_or_default();
    let links = pool.install(|| {
        tasks
            .into_par_iter()
            .map_init(
                || ProteinAligner::with_capacity(max_protein_length, max_protein_length),
                |aligner, task| analyse_gene_task(clusters, task, identity_cutoff, aligner),
            )
            .collect::<Vec<_>>()
    });

    links.into_iter().flatten().collect()
}

#[cfg(not(target_arch = "wasm32"))]
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct AlignmentTask {
    query_cluster: usize,
    target_cluster: usize,
    query_locus: usize,
    query_gene: usize,
}

#[cfg(not(target_arch = "wasm32"))]
fn alignment_tasks(clusters: &[Cluster], first_new_cluster: usize) -> Vec<AlignmentTask> {
    (0..clusters.len())
        .flat_map(|query_cluster_index| {
            let first_target = (query_cluster_index + 1).max(first_new_cluster);
            (first_target..clusters.len()).flat_map(move |target_cluster_index| {
                clusters[query_cluster_index]
                    .loci
                    .iter()
                    .enumerate()
                    .flat_map(move |(query_locus_index, query_locus)| {
                        query_locus
                            .genes
                            .iter()
                            .enumerate()
                            .map(move |(query_gene_index, _)| AlignmentTask {
                                query_cluster: query_cluster_index,
                                target_cluster: target_cluster_index,
                                query_locus: query_locus_index,
                                query_gene: query_gene_index,
                            })
                    })
            })
        })
        .collect()
}

#[cfg(not(target_arch = "wasm32"))]
fn alignment_pool(jobs: usize) -> Result<rayon::ThreadPool, rayon::ThreadPoolBuildError> {
    let threads = if jobs == 0 {
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    } else {
        jobs
    };
    rayon::ThreadPoolBuilder::new().num_threads(threads).build()
}

#[cfg(not(target_arch = "wasm32"))]
fn analyse_gene_task(
    clusters: &[Cluster],
    task: AlignmentTask,
    identity_cutoff: f32,
    aligner: &mut ProteinAligner,
) -> Vec<Link> {
    let query_gene = &clusters[task.query_cluster].loci[task.query_locus].genes[task.query_gene];
    let mut links = Vec::new();

    for (target_locus_index, target_locus) in clusters[task.target_cluster].loci.iter().enumerate()
    {
        for (target_gene_index, target_gene) in target_locus.genes.iter().enumerate() {
            if !can_reach_identity(
                query_gene.translation.len(),
                target_gene.translation.len(),
                identity_cutoff,
            ) {
                continue;
            }
            let ProteinMatch {
                identity,
                similarity,
            } = aligner.compare(
                query_gene.translation.as_bytes(),
                target_gene.translation.as_bytes(),
            );
            if identity >= identity_cutoff {
                links.push(Link {
                    query: GeneRef {
                        cluster: task.query_cluster,
                        locus: task.query_locus,
                        gene: task.query_gene,
                    },
                    target: GeneRef {
                        cluster: task.target_cluster,
                        locus: target_locus_index,
                        gene: target_gene_index,
                    },
                    identity,
                    similarity,
                });
            }
        }
    }

    links
}

/// Compare every gene in `query` with every gene in `target`.
///
/// Returned references are local to this pair: query genes use cluster index
/// zero and target genes use index one. A caller that combines several pair
/// jobs assigns their global cluster indices from its own task metadata.
pub fn analyse_cluster_pair(query: &Cluster, target: &Cluster, identity_cutoff: f32) -> Vec<Link> {
    let mut links = Vec::new();
    let mut aligner =
        ProteinAligner::with_capacity(max_protein_length(query), max_protein_length(target));

    for (query_locus_index, query_locus) in query.loci.iter().enumerate() {
        for (query_gene_index, query_gene) in query_locus.genes.iter().enumerate() {
            for (target_locus_index, target_locus) in target.loci.iter().enumerate() {
                for (target_gene_index, target_gene) in target_locus.genes.iter().enumerate() {
                    let ProteinMatch {
                        identity,
                        similarity,
                    } = if can_reach_identity(
                        query_gene.translation.len(),
                        target_gene.translation.len(),
                        identity_cutoff,
                    ) {
                        aligner.compare(
                            query_gene.translation.as_bytes(),
                            target_gene.translation.as_bytes(),
                        )
                    } else {
                        continue;
                    };

                    if identity >= identity_cutoff {
                        links.push(Link {
                            query: GeneRef {
                                cluster: 0,
                                locus: query_locus_index,
                                gene: query_gene_index,
                            },
                            target: GeneRef {
                                cluster: 1,
                                locus: target_locus_index,
                                gene: target_gene_index,
                            },
                            identity,
                            similarity,
                        });
                    }
                }
            }
        }
    }

    links
}

/// Compare the Cartesian product of two protein slices with the clinker
/// global-alignment metric.
///
/// Callers should keep the slices small enough to form a bounded work tile.
/// The aligner is reused for every pair in the tile.
pub fn analyse_protein_tile(
    query: &[&[u8]],
    target: &[&[u8]],
    identity_cutoff: f32,
) -> Vec<ProteinTileMatch> {
    analyse_protein_tile_with_progress(query, target, identity_cutoff, |_| {})
}

/// As [`analyse_protein_tile`], reporting the number of processed pairs every
/// 100 comparisons and once at completion.
pub fn analyse_protein_tile_with_progress(
    query: &[&[u8]],
    target: &[&[u8]],
    identity_cutoff: f32,
    mut progress: impl FnMut(usize),
) -> Vec<ProteinTileMatch> {
    let mut aligner = ProteinAligner::with_capacity(
        query
            .iter()
            .map(|protein| protein.len())
            .max()
            .unwrap_or_default(),
        target
            .iter()
            .map(|protein| protein.len())
            .max()
            .unwrap_or_default(),
    );
    let mut matches = Vec::new();
    let mut processed = 0;

    for (query_index, query_protein) in query.iter().enumerate() {
        for (target_index, target_protein) in target.iter().enumerate() {
            processed += 1;
            if processed % 100 == 0 {
                progress(processed);
            }
            if !can_reach_identity(query_protein.len(), target_protein.len(), identity_cutoff) {
                continue;
            }
            let ProteinMatch {
                identity,
                similarity,
            } = aligner.compare(query_protein, target_protein);
            if identity >= identity_cutoff {
                matches.push(ProteinTileMatch {
                    query_index,
                    target_index,
                    identity,
                    similarity,
                });
            }
        }
    }

    if processed % 100 != 0 {
        progress(processed);
    }
    matches
}

/// Compare selected pairs from one compact protein collection.
///
/// This is used after an optional candidate filter, so a worker receives each
/// protein once per task instead of expanding selected pairs back into a
/// Cartesian product.
pub fn analyse_protein_pairs(
    proteins: &[&[u8]],
    pairs: &[ProteinPair],
    identity_cutoff: f32,
) -> Vec<ProteinTileMatch> {
    analyse_protein_pairs_with_progress(proteins, pairs, identity_cutoff, |_| {})
}

/// As [`analyse_protein_pairs`], reporting the number of processed pairs
/// every 100 comparisons and once at completion.
pub fn analyse_protein_pairs_with_progress(
    proteins: &[&[u8]],
    pairs: &[ProteinPair],
    identity_cutoff: f32,
    mut progress: impl FnMut(usize),
) -> Vec<ProteinTileMatch> {
    let capacity = proteins
        .iter()
        .map(|protein| protein.len())
        .max()
        .unwrap_or_default();
    let mut aligner = ProteinAligner::with_capacity(capacity, capacity);
    let mut matches = Vec::new();
    let mut processed = 0;

    for &pair in pairs {
        processed += 1;
        if processed % 100 == 0 {
            progress(processed);
        }
        let (Some(query), Some(target)) = (
            proteins.get(pair.query_index),
            proteins.get(pair.target_index),
        ) else {
            continue;
        };
        if !can_reach_identity(query.len(), target.len(), identity_cutoff) {
            continue;
        }
        let ProteinMatch {
            identity,
            similarity,
        } = aligner.compare(query, target);
        if identity >= identity_cutoff {
            matches.push(ProteinTileMatch {
                query_index: pair.query_index,
                target_index: pair.target_index,
                identity,
                similarity,
            });
        }
    }
    if processed % 100 != 0 {
        progress(processed);
    }
    matches
}

/// The largest possible global-alignment identity is the shorter sequence
/// divided by the longer one. Rejecting pairs below this bound is exact.
pub fn can_reach_identity(query_length: usize, target_length: usize, identity_cutoff: f32) -> bool {
    if identity_cutoff <= 0.0 {
        return true;
    }
    let longer = query_length.max(target_length);
    longer != 0 && query_length.min(target_length) as f32 / longer as f32 >= identity_cutoff
}

fn max_protein_length(cluster: &Cluster) -> usize {
    cluster
        .loci
        .iter()
        .flat_map(|locus| &locus.genes)
        .map(|gene| gene.translation.len())
        .max()
        .unwrap_or_default()
}

#[cfg(test)]
mod tests {
    #[cfg(not(target_arch = "wasm32"))]
    use super::{AlignmentTask, alignment_tasks, analyse_genbank_with_options_and_jobs};
    use super::{
        AnalysisOptions, InputFile, ProteinPair, analyse_genbank, analyse_protein_pairs,
        analyse_protein_tile, analyse_protein_tile_with_progress, can_reach_identity,
        parse_input_files_with_progress,
    };
    #[cfg(not(target_arch = "wasm32"))]
    use crate::ParseOptions;
    use crate::{KmerPrefilter, kmer_candidate_pairs, kmer_candidate_pairs_with_progress};

    const FORWARD_CDS: &[u8] =
        br#"LOCUS       FIRST                      9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             1..9
                     /locus_tag="first"
ORIGIN
        1 atggcttaa
//
"#;
    const REVERSE_CDS: &[u8] =
        br#"LOCUS       SECOND                     9 bp    DNA     linear   UNA 01-JAN-2000
FEATURES             Location/Qualifiers
     CDS             complement(1..9)
                     /locus_tag="second"
ORIGIN
        1 ttaagccat
//
"#;

    #[test]
    fn links_similar_genes_from_different_clusters_once() {
        let files = [
            InputFile {
                name: "first.gbk",
                bytes: FORWARD_CDS,
            },
            InputFile {
                name: "second.gbk",
                bytes: REVERSE_CDS,
            },
        ];

        let analysis = analyse_genbank(&files, AnalysisOptions::default()).unwrap();

        assert_eq!(analysis.clusters.len(), 2);
        assert_eq!(analysis.links.len(), 1);
        assert_eq!(analysis.links[0].identity, 1.0);
        assert_eq!(analysis.links[0].query.cluster, 0);
        assert_eq!(analysis.links[0].target.cluster, 1);
        assert_eq!(
            analysis.format_link_summary(),
            "first vs second\n---------------\nQuery  Target  Identity  Similarity\nfirst  second  1.0000    1.0000    "
        );
        assert_eq!(
            analysis.format_link_summary_with_options(2, Some(","), false, false),
            "first,second,1.00,1.00"
        );
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn native_job_counts_preserve_serial_result_order() {
        let files = [
            InputFile {
                name: "first.gbk",
                bytes: FORWARD_CDS,
            },
            InputFile {
                name: "second.gbk",
                bytes: REVERSE_CDS,
            },
            InputFile {
                name: "third.gbk",
                bytes: FORWARD_CDS,
            },
        ];
        let options = AnalysisOptions::default();
        let serial = analyse_genbank(&files, options).unwrap();

        let one_job =
            analyse_genbank_with_options_and_jobs(&files, options, &ParseOptions::default(), 1)
                .unwrap();
        let two_jobs =
            analyse_genbank_with_options_and_jobs(&files, options, &ParseOptions::default(), 2)
                .unwrap();

        assert_eq!(one_job, serial);
        assert_eq!(two_jobs, serial);
    }

    #[cfg(not(target_arch = "wasm32"))]
    #[test]
    fn appended_alignment_tasks_include_only_pairs_with_new_clusters() {
        let files = [
            InputFile {
                name: "first.gbk",
                bytes: FORWARD_CDS,
            },
            InputFile {
                name: "second.gbk",
                bytes: REVERSE_CDS,
            },
            InputFile {
                name: "third.gbk",
                bytes: FORWARD_CDS,
            },
        ];
        let clusters = super::parse_input_files(&files).unwrap();

        assert_eq!(
            alignment_tasks(&clusters, 2),
            vec![
                AlignmentTask {
                    query_cluster: 0,
                    target_cluster: 2,
                    query_locus: 0,
                    query_gene: 0,
                },
                AlignmentTask {
                    query_cluster: 1,
                    target_cluster: 2,
                    query_locus: 0,
                    query_gene: 0,
                },
            ]
        );
    }

    #[test]
    fn parsing_reports_each_completed_file() {
        let files = [
            InputFile {
                name: "first.gbk",
                bytes: FORWARD_CDS,
            },
            InputFile {
                name: "second.gbk",
                bytes: REVERSE_CDS,
            },
        ];
        let mut reports = Vec::new();

        let clusters =
            parse_input_files_with_progress(&files, |processed| reports.push(processed)).unwrap();

        assert_eq!(clusters.len(), 2);
        assert_eq!(reports, vec![1, 2]);
    }

    #[test]
    fn compares_a_bounded_cartesian_product_of_proteins() {
        let query: [&[u8]; 2] = [b"MST", b"AAA"];
        let target: [&[u8]; 2] = [b"GGG", b"MST"];

        let matches = analyse_protein_tile(&query, &target, 0.99);

        assert_eq!(matches.len(), 1);
        assert_eq!(matches[0].query_index, 0);
        assert_eq!(matches[0].target_index, 1);
        assert_eq!(matches[0].identity, 1.0);
    }

    #[test]
    fn length_bound_rejects_pairs_that_cannot_meet_global_identity_cutoff() {
        assert!(can_reach_identity(30, 100, 0.30));
        assert!(!can_reach_identity(29, 100, 0.30));
    }

    #[test]
    fn kmer_prefilter_keeps_only_pairs_with_enough_distinct_shared_words() {
        let query: [&[u8]; 2] = [b"MSTAVK", b"QQQQQQ"];
        let target: [&[u8]; 2] = [b"MSTAVR", b"GGGGGG"];
        let pairs = kmer_candidate_pairs(&query, &target, 0.0, KmerPrefilter::default());

        assert_eq!(pairs.len(), 1);
        assert_eq!(pairs[0].query_index, 0);
        assert_eq!(pairs[0].target_index, 0);
    }

    #[test]
    fn kmer_prefilter_reports_searched_possible_pairs() {
        let query: [&[u8]; 2] = [b"MSTAVK", b"QQQQQQ"];
        let target: [&[u8]; 3] = [b"MSTAVR", b"GGGGGG", b"AAAAAA"];
        let mut reports = Vec::new();

        let _ = kmer_candidate_pairs_with_progress(
            &query,
            &target,
            0.0,
            KmerPrefilter::default(),
            |processed| reports.push(processed),
        );

        assert_eq!(reports, vec![3, 6]);
    }

    #[test]
    fn selected_pair_alignment_does_not_expand_to_a_cartesian_product() {
        let proteins: [&[u8]; 3] = [b"MST", b"AAA", b"GGG"];
        let matches = analyse_protein_pairs(
            &proteins,
            &[ProteinPair {
                query_index: 0,
                target_index: 2,
            }],
            0.99,
        );
        assert!(matches.is_empty());
    }

    #[test]
    fn tile_progress_reports_the_final_processed_pair_count() {
        let proteins: [&[u8]; 2] = [b"MST", b"AAA"];
        let mut reports = Vec::new();
        let _ = analyse_protein_tile_with_progress(&proteins, &proteins, 0.99, |processed| {
            reports.push(processed)
        });
        assert_eq!(reports, vec![4]);
    }
}
