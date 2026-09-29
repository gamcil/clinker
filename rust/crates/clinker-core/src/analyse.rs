use thiserror::Error;

use crate::{Cluster, ProteinMatch, compare_proteins, parse_genbank};

/// An in-memory input boundary shared by native and browser front ends.
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
    /// Format retained links in the same cluster-pair table layout as the
    /// original CLI. This is presentation only: callers that need structured
    /// results should use [`Self::links`] or [`Self::to_plot_data`].
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
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
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

#[derive(Debug, Clone, PartialEq)]
pub struct Analysis {
    pub clusters: Vec<Cluster>,
    pub links: Vec<Link>,
}

/// Parse GenBank inputs and retain cross-cluster gene links at the identity
/// cutoff. Grouping and plot-data serialization follow in later increments.
pub fn analyse_genbank(
    files: &[InputFile<'_>],
    options: AnalysisOptions,
) -> Result<Analysis, AnalysisError> {
    let clusters = files
        .iter()
        .map(|file| {
            parse_genbank(file.name, file.bytes).map_err(|source| AnalysisError::Parse {
                file_name: file.name.to_owned(),
                source,
            })
        })
        .collect::<Result<Vec<_>, _>>()?;

    Ok(Analysis {
        links: cross_cluster_links(&clusters, options.identity_cutoff),
        clusters,
    })
}

fn cross_cluster_links(clusters: &[Cluster], identity_cutoff: f32) -> Vec<Link> {
    let mut links = Vec::new();

    for query_cluster_index in 0..clusters.len() {
        for target_cluster_index in query_cluster_index + 1..clusters.len() {
            for (query_locus_index, query_locus) in
                clusters[query_cluster_index].loci.iter().enumerate()
            {
                for (query_gene_index, query_gene) in query_locus.genes.iter().enumerate() {
                    for (target_locus_index, target_locus) in
                        clusters[target_cluster_index].loci.iter().enumerate()
                    {
                        for (target_gene_index, target_gene) in
                            target_locus.genes.iter().enumerate()
                        {
                            let ProteinMatch {
                                identity,
                                similarity,
                            } = compare_proteins(
                                query_gene.translation.as_bytes(),
                                target_gene.translation.as_bytes(),
                            );

                            if identity >= identity_cutoff {
                                links.push(Link {
                                    query: GeneRef {
                                        cluster: query_cluster_index,
                                        locus: query_locus_index,
                                        gene: query_gene_index,
                                    },
                                    target: GeneRef {
                                        cluster: target_cluster_index,
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
        }
    }

    links
}

#[cfg(test)]
mod tests {
    use super::{AnalysisOptions, InputFile, analyse_genbank};

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
}
