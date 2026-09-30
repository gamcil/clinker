use std::collections::{HashMap, HashSet};

use thiserror::Error;

use crate::align::ProteinAligner;
use crate::{Cluster, ProteinMatch, parse_genbank};

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
                let mut rows = vec![(
                    "Query".to_owned(),
                    "Target".to_owned(),
                    "Identity".to_owned(),
                    "Similarity".to_owned(),
                )];
                rows.extend(links.into_iter().map(|link| {
                    (
                        self.gene_label(link.query).to_owned(),
                        self.gene_label(link.target).to_owned(),
                        format!("{:.4}", link.identity),
                        format!("{:.4}", link.similarity),
                    )
                }));
                let widths = [
                    rows.iter().map(|row| row.0.len()).max().unwrap_or(0),
                    rows.iter().map(|row| row.1.len()).max().unwrap_or(0),
                    rows.iter().map(|row| row.2.len()).max().unwrap_or(0),
                    rows.iter().map(|row| row.3.len()).max().unwrap_or(0),
                ];
                let table = rows
                    .into_iter()
                    .map(|row| {
                        format!(
                            "{:<w0$}  {:<w1$}  {:<w2$}  {:<w3$}",
                            row.0,
                            row.1,
                            row.2,
                            row.3,
                            w0 = widths[0],
                            w1 = widths[1],
                            w2 = widths[2],
                            w3 = widths[3],
                        )
                    })
                    .collect::<Vec<_>>()
                    .join("\n");

                sections.push(format!("{header}\n{}\n{table}", "-".repeat(header.len())));
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
}

/// A stable position within this analysis result.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneRef {
    pub cluster: usize,
    pub locus: usize,
    pub gene: usize,
}

#[derive(Debug, Clone, PartialEq)]
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

/// Conservative settings for the optional, approximate protein k-mer filter.
///
/// Three-residue protein words are conventional for seed-and-extend search,
/// while requiring three distinct hits avoids spending global alignment time on
/// most coincidental single-word matches. These settings can still miss remote
/// homologues, so callers must keep this filter opt-in.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct KmerPrefilter {
    pub kmer_size: usize,
    pub min_shared_kmers: usize,
}

impl Default for KmerPrefilter {
    fn default() -> Self {
        Self {
            kmer_size: 3,
            min_shared_kmers: 3,
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Analysis {
    pub clusters: Vec<Cluster>,
    pub links: Vec<Link>,
}

/// Parse GenBank inputs and retain cross-cluster gene links at the identity
/// cutoff.
pub fn analyse_genbank(
    files: &[InputFile<'_>],
    options: AnalysisOptions,
) -> Result<Analysis, AnalysisError> {
    let clusters = parse_input_files(files)?;

    Ok(Analysis {
        links: cross_cluster_links(&clusters, options.identity_cutoff),
        clusters,
    })
}

/// Parse inputs without comparing proteins.
///
/// Browser callers use this first to obtain cluster metadata, then distribute
/// independent pairwise comparisons to workers.
pub fn parse_input_files(files: &[InputFile<'_>]) -> Result<Vec<Cluster>, AnalysisError> {
    files
        .iter()
        .map(|file| {
            parse_genbank(file.name, file.bytes).map_err(|source| AnalysisError::Parse {
                file_name: file.name.to_owned(),
                source,
            })
        })
        .collect()
}

fn cross_cluster_links(clusters: &[Cluster], identity_cutoff: f32) -> Vec<Link> {
    let mut links = Vec::new();

    for query_cluster_index in 0..clusters.len() {
        for target_cluster_index in query_cluster_index + 1..clusters.len() {
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

/// Return promising pairs from two protein collections using distinct exact
/// peptide k-mers. The returned pairs still require full global alignment.
///
/// This is intentionally a heuristic: unlike the length bound, a k-mer filter
/// can exclude a remote homologue. It is therefore not used by default.
pub fn kmer_candidate_pairs(
    query: &[&[u8]],
    target: &[&[u8]],
    identity_cutoff: f32,
    prefilter: KmerPrefilter,
) -> Vec<ProteinPair> {
    if prefilter.kmer_size == 0 || prefilter.kmer_size > 8 || prefilter.min_shared_kmers == 0 {
        return cartesian_pairs(query, target, identity_cutoff);
    }

    let mut target_index = HashMap::<u64, Vec<usize>>::new();
    for (target_index_value, protein) in target.iter().enumerate() {
        for kmer in distinct_kmers(protein, prefilter.kmer_size) {
            target_index
                .entry(kmer)
                .or_default()
                .push(target_index_value);
        }
    }

    let mut candidates = Vec::new();
    let mut shared_counts = vec![0_usize; target.len()];
    let mut touched_targets = Vec::new();
    for (query_index, protein) in query.iter().enumerate() {
        for kmer in distinct_kmers(protein, prefilter.kmer_size) {
            if let Some(targets) = target_index.get(&kmer) {
                for &target_index in targets {
                    if shared_counts[target_index] == 0 {
                        touched_targets.push(target_index);
                    }
                    shared_counts[target_index] += 1;
                }
            }
        }
        for target_index in touched_targets.drain(..) {
            let count = std::mem::take(&mut shared_counts[target_index]);
            if count >= prefilter.min_shared_kmers
                && can_reach_identity(protein.len(), target[target_index].len(), identity_cutoff)
            {
                candidates.push(ProteinPair {
                    query_index,
                    target_index,
                });
            }
        }
    }
    candidates.sort_unstable_by_key(|pair| (pair.query_index, pair.target_index));
    candidates
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

fn cartesian_pairs(query: &[&[u8]], target: &[&[u8]], identity_cutoff: f32) -> Vec<ProteinPair> {
    query
        .iter()
        .enumerate()
        .flat_map(|(query_index, query_protein)| {
            target
                .iter()
                .enumerate()
                .filter_map(move |(target_index, target_protein)| {
                    can_reach_identity(query_protein.len(), target_protein.len(), identity_cutoff)
                        .then_some(ProteinPair {
                            query_index,
                            target_index,
                        })
                })
        })
        .collect()
}

fn distinct_kmers(protein: &[u8], kmer_size: usize) -> HashSet<u64> {
    protein
        .windows(kmer_size)
        .map(|kmer| {
            kmer.iter()
                .fold(0_u64, |key, &amino_acid| (key << 8) | u64::from(amino_acid))
        })
        .collect()
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
    use super::{
        AnalysisOptions, InputFile, KmerPrefilter, ProteinPair, analyse_genbank,
        analyse_protein_pairs, analyse_protein_tile, analyse_protein_tile_with_progress,
        can_reach_identity, kmer_candidate_pairs,
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
