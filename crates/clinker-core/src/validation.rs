use std::collections::HashSet;

use thiserror::Error;

use crate::{Analysis, GeneRef};

#[derive(Debug, Error)]
pub enum ValidationError {
    #[error("cluster {0} contains no loci")]
    EmptyCluster(usize),
    #[error("locus ({cluster}, {locus}) must have start < end")]
    LocusCoordinates { cluster: usize, locus: usize },
    #[error("gene {0:?} must have start < end and be contained within its locus")]
    GeneCoordinates(GeneRef),
    #[error("gene {0:?} must have strand -1 or 1")]
    GeneStrand(GeneRef),
    #[error("link {link} has an invalid {endpoint} reference {reference:?}")]
    GeneReference {
        link: usize,
        endpoint: &'static str,
        reference: GeneRef,
    },
    #[error("link {link}: {reason}")]
    Link { link: usize, reason: &'static str },
    #[error("cluster order must include every cluster exactly once")]
    ClusterOrder,
}

impl Analysis {
    pub fn validate(&self) -> Result<(), ValidationError> {
        for (cluster_index, cluster) in self.clusters.iter().enumerate() {
            if cluster.loci.is_empty() {
                return Err(ValidationError::EmptyCluster(cluster_index));
            }
            for (locus_index, locus) in cluster.loci.iter().enumerate() {
                if locus.start >= locus.end {
                    return Err(ValidationError::LocusCoordinates {
                        cluster: cluster_index,
                        locus: locus_index,
                    });
                }
                for (gene_index, gene) in locus.genes.iter().enumerate() {
                    let reference = GeneRef {
                        cluster: cluster_index,
                        locus: locus_index,
                        gene: gene_index,
                    };
                    if gene.start >= gene.end || gene.start < locus.start || gene.end > locus.end {
                        return Err(ValidationError::GeneCoordinates(reference));
                    }
                    if !matches!(gene.strand, -1 | 1) {
                        return Err(ValidationError::GeneStrand(reference));
                    }
                }
            }
        }
        let mut pairs = HashSet::new();
        for (index, link) in self.links.iter().enumerate() {
            for (endpoint, reference) in [("query", link.query), ("target", link.target)] {
                if self
                    .clusters
                    .get(reference.cluster)
                    .and_then(|cluster| cluster.loci.get(reference.locus))
                    .and_then(|locus| locus.genes.get(reference.gene))
                    .is_none()
                {
                    return Err(ValidationError::GeneReference {
                        link: index,
                        endpoint,
                        reference,
                    });
                }
            }
            let reason = if link.query.cluster >= link.target.cluster {
                Some("query cluster must precede target cluster")
            } else if !(0.0..=1.0).contains(&link.identity)
                || !(link.identity..=1.0).contains(&link.similarity)
            {
                Some("scores must satisfy 0 <= identity <= similarity <= 1")
            } else if !pairs.insert((link.query, link.target)) {
                Some("duplicate gene pair")
            } else {
                None
            };
            if let Some(reason) = reason {
                return Err(ValidationError::Link {
                    link: index,
                    reason,
                });
            }
        }
        Ok(())
    }

    pub fn validate_order(&self, order: &[usize]) -> Result<(), ValidationError> {
        if order.len() != self.clusters.len() {
            return Err(ValidationError::ClusterOrder);
        }
        let mut seen = vec![false; self.clusters.len()];
        for &cluster in order {
            if cluster >= seen.len() || std::mem::replace(&mut seen[cluster], true) {
                return Err(ValidationError::ClusterOrder);
            }
        }
        Ok(())
    }
}
