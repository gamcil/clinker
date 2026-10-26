use std::collections::HashSet;

use serde::{Deserialize, Serialize};
use thiserror::Error;

#[cfg(not(target_arch = "wasm32"))]
use rayon::prelude::*;

use crate::{
    AlignmentCounts, Cluster, KmerPrefilter, ParseOptions, PreparedAlignment, ProteinBatchAligner,
    parse_genbank_with_options, parse_gff3,
};

#[derive(Debug, Clone, Copy)]
pub struct InputFile<'a> {
    pub name: &'a str,
    pub bytes: &'a [u8],
}

#[derive(Debug, Clone, Copy)]
pub struct AnalysisOptions {
    pub identity_cutoff: f32,
    /// Optional approximate candidate selection before alignment.
    pub prefilter: Option<KmerPrefilter>,
}

/// Work performed by the alignment workflow, after length and optional k-mer filtering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct AlignmentStats {
    pub aligned_pairs: usize,
    pub counts: AlignmentCounts,
}

impl Default for AnalysisOptions {
    fn default() -> Self {
        Self {
            identity_cutoff: 0.30,
            prefilter: None,
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
#[derive(Debug, Clone, Copy, PartialEq, Eq, Deserialize, Serialize)]
#[serde(rename_all = "camelCase")]
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
            AnalysisOptions {
                identity_cutoff,
                ..AnalysisOptions::default()
            },
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
        self.append_clusters_with_options_and_jobs(
            clusters,
            AnalysisOptions {
                identity_cutoff,
                ..AnalysisOptions::default()
            },
            jobs,
        )
        .map(|_| ())
    }

    /// Append clusters using the same native scheduler with or without prefiltering.
    #[cfg(not(target_arch = "wasm32"))]
    pub fn append_clusters_with_options_and_jobs(
        &mut self,
        clusters: Vec<Cluster>,
        options: AnalysisOptions,
        jobs: usize,
    ) -> Result<AlignmentStats, AnalysisError> {
        let pool = alignment_pool(jobs)?;
        let first_new_cluster = self.clusters.len();
        self.clusters.extend(clusters);
        let (links, stats) =
            cross_cluster_links_parallel_from(&self.clusters, first_new_cluster, options, &pool);
        self.links.extend(links);
        Ok(stats)
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
        links: cross_cluster_links_from(&clusters, 0, options),
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
        links: cross_cluster_links_parallel_from(&clusters, 0, options, &pool).0,
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

fn collect_batch_links(
    results: impl IntoIterator<Item = (Vec<Link>, usize)>,
    counts: AlignmentCounts,
) -> (Vec<Link>, AlignmentStats) {
    let mut links = Vec::new();
    let mut stats = AlignmentStats {
        aligned_pairs: 0,
        counts,
    };
    for (batch_links, aligned_pairs) in results {
        links.extend(batch_links);
        stats.aligned_pairs += aligned_pairs;
    }
    // Tile boundaries must not change link order or stable IDs in native output.
    links.sort_unstable_by_key(|link| {
        (
            link.query.cluster,
            link.target.cluster,
            link.query.locus,
            link.query.gene,
            link.target.locus,
            link.target.gene,
        )
    });
    (links, stats)
}

fn cross_cluster_links_from(
    clusters: &[Cluster],
    first_new_cluster: usize,
    options: AnalysisOptions,
) -> Vec<Link> {
    let inputs = PreparedAlignment::new(clusters, first_new_cluster, options, |_, _| {});
    let mut aligner = ProteinBatchAligner::default();
    collect_batch_links(
        inputs
            .batches
            .iter()
            .map(|batch| inputs.run(batch, &mut aligner, options.identity_cutoff)),
        inputs.counts,
    )
    .0
}

#[cfg(not(target_arch = "wasm32"))]
fn cross_cluster_links_parallel_from(
    clusters: &[Cluster],
    first_new_cluster: usize,
    options: AnalysisOptions,
    pool: &rayon::ThreadPool,
) -> (Vec<Link>, AlignmentStats) {
    let inputs = PreparedAlignment::new(clusters, first_new_cluster, options, |_, _| {});
    let results = pool.install(|| {
        inputs
            .batches
            .par_iter()
            .map_init(ProteinBatchAligner::default, |aligner, batch| {
                inputs.run(batch, aligner, options.identity_cutoff)
            })
            .collect::<Vec<_>>()
    });
    collect_batch_links(results, inputs.counts)
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

/// The largest possible global-alignment identity is the shorter sequence
/// divided by the longer one. Rejecting pairs below this bound is exact.
pub fn can_reach_identity(query_length: usize, target_length: usize, identity_cutoff: f32) -> bool {
    if identity_cutoff <= 0.0 {
        return true;
    }
    let longer = query_length.max(target_length);
    longer != 0 && query_length.min(target_length) as f32 / longer as f32 >= identity_cutoff
}

#[cfg(test)]
mod tests {
    #[cfg(not(target_arch = "wasm32"))]
    use super::analyse_genbank_with_options_and_jobs;
    use super::{
        AnalysisOptions, InputFile, ProteinPair, analyse_genbank, can_reach_identity,
        parse_input_files_with_progress,
    };
    #[cfg(not(target_arch = "wasm32"))]
    use crate::ParseOptions;
    use crate::{
        KmerPrefilter, ProteinBatch, kmer_candidate_pairs, kmer_candidate_pairs_with_progress,
    };

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
    fn prefiltered_workflow_preserves_references_order_and_incremental_links() {
        let template = super::parse_input_files(&[InputFile {
            name: "first.gbk",
            bytes: FORWARD_CDS,
        }])
        .unwrap()
        .remove(0);
        let mut clusters = vec![template; 3];
        for cluster in &mut clusters {
            cluster.loci.push(cluster.loci[0].clone());
            cluster.loci[0].genes[0].translation = "MSTAVKLMDEFGH".into();
            cluster.loci[1].genes[0].translation = "MSTAVRLMDEFGH".into();
            let mut short = cluster.loci[1].genes[0].clone();
            short.translation = "MS".into();
            cluster.loci[1].genes.push(short);
        }
        for (first_new, expected) in [(0, 27), (2, 18), (3, 0)] {
            let prepared = super::PreparedAlignment::new(
                &clusters,
                first_new,
                AnalysisOptions::default(),
                |_, _| {},
            );
            assert_eq!(prepared.counts.possible_pairs, expected);
            assert_eq!(prepared.counts.planned_pairs, expected);
            assert_eq!(prepared.counts.batch_count, prepared.batches.len());
        }
        let unfiltered = super::cross_cluster_links_from(&clusters, 0, AnalysisOptions::default());
        for prefilter in [
            KmerPrefilter::default(),
            KmerPrefilter {
                kmer_size: 5,
                min_shared_kmers: 8,
            },
            KmerPrefilter {
                kmer_size: 8,
                min_shared_kmers: 100,
            },
        ] {
            let options = AnalysisOptions {
                prefilter: Some(prefilter),
                ..AnalysisOptions::default()
            };
            let prepared = super::PreparedAlignment::new(&clusters, 0, options, |_, _| {});
            assert_eq!(
                prepared.counts.planned_pairs,
                prepared.candidates.as_ref().unwrap().len()
            );
            assert_eq!(
                prepared.counts.planned_pairs,
                prepared
                    .batches
                    .iter()
                    .map(ProteinBatch::pair_count)
                    .sum::<usize>()
            );
            let expected = unfiltered
                .iter()
                .filter(|link| {
                    let protein = |r: super::GeneRef| {
                        clusters[r.cluster].loci[r.locus].genes[r.gene]
                            .translation
                            .as_bytes()
                    };
                    !kmer_candidate_pairs(
                        &[protein(link.query)],
                        &[protein(link.target)],
                        options.identity_cutoff,
                        prefilter,
                    )
                    .is_empty()
                })
                .cloned()
                .collect::<Vec<_>>();
            assert_eq!(
                super::cross_cluster_links_from(&clusters, 0, options),
                expected
            );
            for jobs in [1, 2] {
                let mut analysis = super::Analysis {
                    clusters: Vec::new(),
                    links: Vec::new(),
                };
                let stats = analysis
                    .append_clusters_with_options_and_jobs(clusters.clone(), options, jobs)
                    .unwrap();
                assert_eq!(analysis.links, expected);
                assert!(stats.aligned_pairs >= analysis.links.len());
                if prefilter.min_shared_kmers == 100 {
                    assert_eq!(stats.aligned_pairs, 0);
                }
                let mut incremental = super::Analysis {
                    clusters: Vec::new(),
                    links: Vec::new(),
                };
                incremental
                    .append_clusters_with_options_and_jobs(clusters[..2].to_vec(), options, jobs)
                    .unwrap();
                let existing = incremental.links.clone();
                incremental
                    .append_clusters_with_options_and_jobs(clusters[2..].to_vec(), options, jobs)
                    .unwrap();
                assert_eq!(incremental.links, expected);
                assert!(incremental.links.starts_with(&existing));
            }
        }
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

        let batch = ProteinBatch::Cartesian {
            query: 0..query.len(),
            target: 0..target.len(),
        };
        let matches = super::ProteinBatchAligner::default()
            .analyse_pairs(&query, &target, batch.pairs(&[]), 0.99, |_| {})
            .matches;

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
        let matches = super::ProteinBatchAligner::default()
            .analyse_pairs(
                &proteins,
                &proteins,
                [ProteinPair {
                    query_index: 0,
                    target_index: 2,
                }],
                0.99,
                |_| {},
            )
            .matches;
        assert!(matches.is_empty());
    }

    #[test]
    fn tile_progress_reports_the_final_processed_pair_count() {
        let proteins: [&[u8]; 2] = [b"MST", b"AAA"];
        let mut reports = Vec::new();
        let batch = ProteinBatch::Cartesian {
            query: 0..proteins.len(),
            target: 0..proteins.len(),
        };
        let _ = super::ProteinBatchAligner::default().analyse_pairs(
            &proteins,
            &proteins,
            batch.pairs(&[]),
            0.99,
            |processed| reports.push(processed),
        );
        assert_eq!(reports, vec![4]);
    }
}
