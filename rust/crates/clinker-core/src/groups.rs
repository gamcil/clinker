//! Homology groups built from retained gene-to-gene links.

use std::collections::BTreeMap;

use crate::{GeneRef, Link};

/// One connected component of genes joined by retained homology links.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GeneGroup {
    pub genes: Vec<GeneRef>,
}

/// Build homology groups from gene-gene links.
/// A group contains every gene reachable through one or more links.
/// Genes that have no retained link are deliberately omitted.
pub fn build_gene_groups(links: &[Link]) -> Vec<GeneGroup> {
    let mut gene_indices = BTreeMap::new();
    let mut genes = Vec::new();
    for link in links {
        for gene in [link.query, link.target] {
            if !gene_indices.contains_key(&gene) {
                let index = genes.len();
                gene_indices.insert(gene, index);
                genes.push(gene);
            }
        }
    }
    let mut sets = DisjointSet::new(genes.len());
    for link in links {
        sets.union(gene_indices[&link.query], gene_indices[&link.target]);
    }
    let mut components = BTreeMap::<usize, Vec<GeneRef>>::new();
    for (index, gene) in genes.into_iter().enumerate() {
        components.entry(sets.find(index)).or_default().push(gene);
    }
    components
        .into_values()
        .map(|genes| GeneGroup { genes })
        .collect()
}

#[derive(Debug)]
struct DisjointSet {
    parent: Vec<usize>,
    rank: Vec<u8>,
}

impl DisjointSet {
    fn new(len: usize) -> Self {
        Self {
            parent: (0..len).collect(),
            rank: vec![0; len],
        }
    }

    fn find(&mut self, index: usize) -> usize {
        if self.parent[index] != index {
            self.parent[index] = self.find(self.parent[index]);
        }
        self.parent[index]
    }

    fn union(&mut self, one: usize, two: usize) {
        let (one_root, two_root) = (self.find(one), self.find(two));
        if one_root == two_root {
            return;
        }
        if self.rank[one_root] < self.rank[two_root] {
            self.parent[one_root] = two_root;
        } else {
            self.parent[two_root] = one_root;
            if self.rank[one_root] == self.rank[two_root] {
                self.rank[one_root] += 1;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::build_gene_groups;
    use crate::{GeneRef, Link};

    fn gene(cluster: usize, gene: usize) -> GeneRef {
        GeneRef {
            cluster,
            locus: 0,
            gene,
        }
    }

    fn link(query: GeneRef, target: GeneRef) -> Link {
        Link {
            query,
            target,
            identity: 1.0,
            similarity: 1.0,
        }
    }

    #[test]
    fn joins_transitive_links_and_keeps_disconnected_components_separate() {
        let groups = build_gene_groups(&[
            link(gene(0, 0), gene(1, 0)),
            link(gene(1, 0), gene(2, 0)),
            link(gene(0, 1), gene(1, 1)),
        ]);

        assert_eq!(groups.len(), 2);
        assert_eq!(groups[0].genes, vec![gene(0, 0), gene(1, 0), gene(2, 0)]);
        assert_eq!(groups[1].genes, vec![gene(0, 1), gene(1, 1)]);
    }
}
