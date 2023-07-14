#!/usr/bin/env python3

"""
translation protein sequences in Clusters or lists of Proteins.

Cameron Gilchrist
"""

import logging
import uuid
import csv

from collections import defaultdict, OrderedDict
from functools import partial
from itertools import combinations, product
from multiprocessing import Pool
from typing import Dict, List, Optional

import numpy as np

from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform

from Bio import Align, __version__ as biopython_version
from Bio.Align import substitution_matrices

from disjoint_set import DisjointSet

from clinker.formatters import format_alignment, format_globaligner
from clinker.classes import Serializer, Cluster, Locus, Gene, load_child, load_children


LOG = logging.getLogger(__name__)


def align_clusters(*args, cutoff=0.3, aligner_config=None, jobs=None):
    """Convenience function for directly aligning Cluster object/s.

    Initialises a Globaligner, adds Cluster/s, then runs alignments
    and returns the Globaligner.

    Args:
        *args: Cluster or list of Protein objects
        aligner_config (dict): keyword arguments to use when setting
                               up the BioPython.PairwiseAligner object
        cutoff (float): decimal identity cutoff for saving an alignment
        jobs (int, optional): number of jobs to run alignment in parallel
    Returns:
        aligner (Globaligner): instance of Globaligner class which
                                  contains all cluster alignments
    e.g.
        align_sequence_groups(cluster1, cluster2, ..., clusterN)
    """
    aligner = Globaligner(aligner_config)
    aligner.add_clusters(*args)
    if len(args) == 1:
        LOG.info("Only one cluster given, skipping alignment")
    else:
        aligner.align_stored_clusters(cutoff, jobs=jobs)
    return aligner


def assign_groups(links, threshold=0.3):
    """Groups sequences in alignment links by single-linkage."""
    groups = []
    for link in links:
        if link.identity < threshold:
            continue
        found = False
        for (i, group) in enumerate(groups):
            if link.query in group or link.target in group:
                found = True
            if found:
                for gene in [link.query, link.target]:
                    if gene not in group:
                        group.append(gene)
                        gene._group = i
                break
        if not found:
            groups.append([link.query, link.target])
            index = len(groups) - 1
            link.query._group = index
            link.target._group = index


def get_pairs(cluster):
    """Gets all contiguous pairs of homology groups in a cluster."""
    pairs = []
    for locus in cluster.loci:
        total = len(locus.genes) - 1
        pairs.extend(
            (gene._group for gene in locus.genes[i:i+2])
            for i in range(total)
        )
    return pairs


def compare_pairs(one, two):
    """Compares two collections of contiguous group pairs.

    Gets common elements (i.e. intersection) between each list, and then
    finds the minimum number of occurrences of the elements in either,
    such that shared duplicate pairs will be included in the total.
    """
    total = 0
    for pair in set(one).intersection(two):
        total += min(one.count(pair), two.count(pair))
    return total


def compute_identity(alignment):
    """Calculates sequence identity/similarity of a BioPython alignment object."""

    similar_acids = [
        {"G", "A", "V", "L", "I"},
        {"F", "Y", "W"},
        {"C", "M"},
        {"S", "T"},
        {"K", "R", "H"},
        {"D", "E", "N", "Q"},
        {"P"},
    ]
                
    if biopython_version >= "1.80":
        # Default format changed as of BioPython v1.80
        # https://github.com/biopython/biopython/issues/4183
        one, two = alignment
    else:
        one, _, two, _ = str(alignment).split("\n")

    length = len(one)

    # Amino acid similarity groups
    matches, similar = 0, 0
    for i in range(length):
        if one[i] == two[i]:
            # Check for gap columns
            if one[i] not in {"-", "."}:
                matches += 1
            else:
                length -= 1
        else:
            # If not identical, check if similar
            for group in similar_acids:
                if one[i] in group and two[i] in group:
                    similar += 1
                    break

    # identity = matches / length - gaps
    # similarity = (matches + similarities) / length - gaps
    return matches / length, (matches + similar) / length


def extend_matrix_alphabet(matrix, codes='BXZJUO'):
    """Extends the alphabet of a given substitution matrix.

    Primarily for adding extended IUPAC codes to a matrix which does
    not contain them (e.g. BLOSUM62), resulting in a ValueError
    being thrown during sequence alignment.
    """
    missing_codes = set(codes).difference(matrix.alphabet)
    if missing_codes:
        missing_codes = ''.join(missing_codes)
        matrix = matrix.select(matrix.alphabet + missing_codes)
    return matrix


class Globaligner(Serializer):
    """Performs and stores alignments.

    Parameters:
        aligner (Bio.Align.PairwiseAligner): Sqeuence aligner
        alignments (list): Alignments generated by Globaligner
        clusters (dict): Ordered dictionary of Clusters keyed on name
        _alignment_indices (dict): indices of Alignments in _alignments
            stored using Cluster.name attributes as keys
        _cluster_names (dict): tuples of Cluster.name attributes stored using
            _alignment indices as keys
    """

    aligner_default = {
        "mode": "global",
        "substitution_matrix": "BLOSUM62",
        "open_gap_score": -10,
        "extend_gap_score": -0.5,
    }

    def __init__(self, aligner_config=None):
        # Lookup dictionaries
        self._genes = {}
        self._loci = {}
        self._links = {}
        self._alignment_indices = defaultdict(dict)
        self._cluster_names = defaultdict(dict)

        self.alignments = {}
        self.groups = []
        self.clusters = OrderedDict()

        if aligner_config is None:
            self.aligner_config = self.aligner_default.copy()
        else:
            self.aligner_config = aligner_config

    def to_dict(self):
        """Serialises the Globaligner instance to dict.

        Schema:
            {
                order: [],
                clusters: {},
                loci: {},
                genes: {},
                alignments: [],
                links: {},
            }

        Where each child dictionary holds serialised Python objects keyed
        by their UIDs. When from_dict() is called, they are used to store
        real references between objects (e.g. Link query/target attributes
        are Gene objects).
        """
        return {
            "order": list(self.clusters),
            "clusters": {
                uid: cluster.to_dict(uids_only=True)
                for uid, cluster in self.clusters.items()
            },
            "loci": {
                uid: locus.to_dict(uids_only=True)
                for uid, locus in self._loci.items()
            },
            "genes": {
                uid: gene.to_dict()
                for uid, gene in self._genes.items()
            },
            "groups": [group.to_dict() for group in self.groups],
            "alignments": {
                uid: alignment.to_dict(uids_only=True)
                for uid, alignment in self.alignments.items()
            },
            "links": {
                uid: link.to_dict()
                for uid, link in self._links.items()
            },
        }

    @classmethod
    def from_dict(cls, d):
        """Loads a Globaligner instance from dict generated by to_dict().

        First, loads all clinker objects back into memory (cluster, locus, gene)
        and restores their hierarchical structure. Alignments are restored in the same way.
        Finally, rebuilds lookup dictionaries used by the Globaligner class.
        """
        ga = Globaligner()

        # Reconstruct cluster -> locus -> gene hierarchy
        for cluster_uid in d["order"]:
            cluster = Cluster.from_dict(d["clusters"][cluster_uid])
            for locus_idx, locus_uid in enumerate(cluster.loci):
                locus = Locus.from_dict(d["loci"][locus_uid])
                for gene_idx, gene_uid in enumerate(locus.genes):
                    gene = Gene.from_dict(d["genes"][gene_uid])
                    locus.genes[gene_idx] = gene
                    ga._genes[gene_uid] = gene
                cluster.loci[locus_idx] = locus
                ga._loci[locus_uid] = locus
            ga.clusters[cluster_uid] = cluster

        # Reconstruct cluster alignments
        for alignment_uid, alignment in d["alignments"].items():
            aln = Alignment.from_dict(alignment)
            aln.query = ga.clusters[aln.query]
            aln.target = ga.clusters[aln.target]

            # Form Link objects and add to lookup dict
            for idx, uid in enumerate(aln.links):
                link = Link.from_dict(d["links"][uid])
                link.query = ga._genes[link.query.uid]
                link.target = ga._genes[link.target.uid]
                aln.links[idx] = link
                ga._links[uid] = link

            # Update Alignment object lookup dictionaries
            ga.alignments[aln.uid] = aln
            ga._alignment_indices[aln.query.uid][aln.target.uid] = aln.uid
            ga._alignment_indices[aln.target.uid][aln.query.uid] = aln.uid
            ga._cluster_names[aln.uid] = (aln.query.uid, aln.target.uid)

        for group in d["groups"]:
            gr = Group.from_dict(group)
            ga.groups.append(gr)

        return ga

    def __str__(self):
        """Print all alignments currently stored in the instance."""
        return self.format()

    def format(
        self,
        delimiter=None,
        decimals=4,
        alignment_headers=True,
        link_headers=False,
    ):
        return format_globaligner(
            self,
            decimals=decimals,
            delimiter=delimiter,
            alignment_headers=alignment_headers,
            link_headers=link_headers,
        )

    def to_data(self, i=0.5, method="ward", use_file_order=False):
        """Formats Globaligner as plottable data set.

        Assign unique indices to all clusters, loci, genes
        """
        clusters = [cluster.to_dict() for cluster in self.clusters.values()]
        return {
            "clusters": clusters if use_file_order else [
                clusters[i] for i in self.order(i=i, method=method)
            ],
            "links": [
                link.to_dict()
                for alignment in self.alignments.values()
                for link in alignment.links
            ],
            "groups": [group.to_dict() for group in self.groups]
        }

    @property
    def gene_labels(self):
        labels = set()
        for cluster in self.clusters.values():
            for locus in cluster.loci:
                for gene in locus.genes:
                    labels.update(gene.names)
        return labels

    def add_clusters(self, *clusters):
        """Adds new Cluster object/s to the Globaligner.

        Parameters:
            clusters (list): variable number of Cluster objects
        """
        for cluster in clusters:
            if not isinstance(cluster, Cluster):
                raise NotImplementedError("Expected Cluster object")
            self.clusters[cluster.uid] = cluster
            for locus in cluster.loci:
                self._loci[locus.uid] = locus
                for gene in locus.genes:
                    self._genes[gene.uid] = gene

    @staticmethod
    def _align_clusters(config, one, two, cutoff=0.3):
        """Constructs a cluster alignment using the given configuration."""
        LOG.info("%s vs %s", one.name, two.name)

        aligner = Align.PairwiseAligner()

        # Select the substitution matrix.
        # Defaults to BLOSUM62 when none or invalid matrix specified.
        matrix = config.pop("substitution_matrix", "BLOSUM62")
        if matrix not in substitution_matrices.load():
            LOG.warning("Invalid substitution matrix '(%s)', defaulting to BLOSUM62", matrix)
            matrix = "BLOSUM62"
        aligner.substitution_matrix = substitution_matrices.load(matrix)

        # ValueError is thrown during sequence alignment when a letter
        # in the sequence is not found in the substitution matrix.
        # Extended IUPAC codes (BXZJUO) are added to mitigate this.
        aligner.substitution_matrix = extend_matrix_alphabet(
            aligner.substitution_matrix,
            codes='BXZJUO-.',
        )

        for k, v in config.items():
            setattr(aligner, k, v)

        alignment = Alignment(query=one, target=two)
        for locusA, locusB in product(one.loci, two.loci):
            for geneA, geneB in product(locusA.genes, locusB.genes):
                if not geneA.translation or not geneB.translation:
                    continue
                try:
                    aln = aligner.align(geneA.translation.strip(),
                                        geneB.translation.strip())
                except:
                    raise
                identity, similarity = compute_identity(aln[0])
                if identity < cutoff:
                    continue
                alignment.add_link(geneA, geneB, identity, similarity)
        return alignment

    def align_clusters(self, one, two, cutoff=0.3):
        """Constructs a cluster alignment using aligner config in the Globaligner."""
        return self._align_clusters(self.aligner_config, one, two, cutoff=cutoff)

    def align_stored_clusters(self, cutoff=0.3, jobs=None):
        """Aligns clusters stored in the Globaligner."""

        pairs_to_align = []
        for one, two in combinations(self.clusters.values(), 2):
            if self._alignment_indices[one.uid].get(two.uid):
                LOG.debug("Skipping %s vs %s", one.name, two.name)
            else:
                pairs_to_align.append((one, two))

        with Pool(jobs) as pool:
            _align_clusters = partial(
                self._align_clusters,
                self.aligner_config,
                cutoff=cutoff
            )
            alignments = pool.starmap(_align_clusters, pairs_to_align)

        for alignment in alignments:
            self.add_alignment(alignment)

    def get_gene_uid(self, gene: str) -> Optional[str]:
        """Find UID of a gene stored in the Globaligner given some label."""
        for uid, _gene in self._genes.items():
            if any(gene == value for value in _gene.names.values()):
                return uid

    def get_gene_uids(self, genes: List[str]) -> Optional[List[str]]:
        """Find UIDs of a collection of genes stored in the Globaligner."""
        uids = []
        for gene in genes:
            uid = self.get_gene_uid(gene)
            if not uid:
                LOG.warning("Failed to find gene: %s", gene)
                continue
            uids.append(uid)
        return uids

    def build_gene_groups(
            self,
            functions: Optional[Dict[str, List[str]]]=None,
            colours: Optional[Dict[str, str]]=None
        ) -> None:
        """Builds gene groups based on functions and stored gene-gene links.

        `functions` maps genes to user-assigned functions; keys should correspond
        to the `label` property of Gene objects. If specified, groups will first
        be built from these functional groups and extended using gene-gene links.
        If not, all groups will be built purely from gene-gene links.

        e.g. If genes A and B are both assigned as category Z, and have links to
             genes C and D, respectively, the result will be a single group Z,
             containing genes A, B, C and D. Group Z is first made with genes A and B,
             then extended to include C (because of link to A) and D (link to B).
        """
        self.groups = []
        if functions:
            for function, genes in functions.items():
                uids = self.get_gene_uids(genes)
                group = Group(label=function, genes=set(uids))
                if colours and function in colours:
                    group.colour = colours[function]
                self.groups.append(group)
        if not self._links:
            for group in self.groups:
                group.genes = list(group.genes)
            return
        ds = DisjointSet()
        for link in self._links.values():
            ds.union(link.query.uid, link.target.uid)
        for genes in ds.itersets():
            merged = False
            for group in self.groups:
                if not genes.isdisjoint(group.genes):
                    group.genes.update(genes)
                    merged = True
                    break
            if not merged:
                group = Group(label=f"Group {len(self.groups)}", genes=set(genes))
                self.groups.append(group)
        for group in self.groups:
            group.genes = list(group.genes)

    def configure_aligner(self, **kwargs):
        """Change properties on the BioPython.PairwiseAligner object.

        This function takes any keyword argument and assumes
        they correspond to valid properties on the PairwiseAligner.
        Refer to BioPython documentation for these.
        """
        valid_attributes = {
            x for x in dir(Align.PairwiseAligner)
            if not x.startswith("_")
        }
        invalid_keys = set(kwargs.keys()) - valid_attributes
        if invalid_keys:
            raise ValueError(
                f'invalid attributes for the BioPython Align.PairwiseAligner'
                f"class: { ', '.join(map(repr, sorted(invalid_keys))) }"
            )
        self.aligner_config.update(kwargs)

    def add_alignment(self, alignment, overwrite=False):
        """Adds a new cluster alignment to the Globaligner.

        self._alignment_indices allows for Alignment indices to be
        retrieved from cluster names, regardless of order.

        self._cluster_names allows for Cluster names to be retrieved
        given the index of an Alignment in self.alignments
        """
        # Save Cluster object if not already stored
        q = alignment.query
        t = alignment.target
        self.add_clusters(q, t)

        # Overwrite previous alignment between these clusters if one exists
        previous = self._alignment_indices[q.uid].get(t.uid)

        # Clear up any old links, add the new ones
        if previous:
            previous = self.alignments.pop(previous)
            for link in previous.links:
                self._links.pop(link.uid)
        for link in alignment.links:
            self._links[link.uid] = link

        # Update mapping dictionaries and save Alignment
        self.alignments[alignment.uid] = alignment
        self._alignment_indices[q.uid][t.uid] = alignment.uid
        self._alignment_indices[t.uid][q.uid] = alignment.uid
        self._cluster_names[alignment.uid] = (q.uid, t.uid)

    def get_alignment(self, one, two):
        """Retrieves an Alignment corresponding to two Cluster objects.

        Parameters:
            one (str): Name of first cluster
            two (str): Name of second cluster
        Returns:
            Alignment object for the specified Clusters
        """
        uid = self._alignment_indices[one][two]
        return self.alignments[uid]

    def synteny(self, one, two, i=0.5):
        """Calculates a synteny score between two clusters.

        Based on antiSMASH/MultiGeneBlast implementation:
            S = h + i*s
        where:
            h = #homologues over minimum identity/coverage threshold
            s = #contiguous gene pairs
            i = weighting factor for s

        Except instead of counting number of homologues, we use a cumulative
        identity value of homologues in each cluster.
        """
        alignment = self.get_alignment(one, two)
        if not alignment.links:
            return 0
        homology = sum(link.identity for link in alignment.links)
        assign_groups(alignment.links)
        one_cluster = self.clusters[one]
        two_cluster = self.clusters[two]
        one_pairs = get_pairs(one_cluster)
        two_pairs = get_pairs(two_cluster)
        contiguity = compare_pairs(one_pairs, two_pairs)
        return homology + i * contiguity

    def matrix(self, i=0.5, normalise=False, as_distance=False):
        """Generates a synteny score matrix of all aligned clusters.

        Arguments:
            i (float): Weighting of gene pair contiguity in synteny scores
            normalise (bool): Normalise the matrix (i.e. 0 to 1)
            as_distance (bool): Convert to distance matrix
        Returns:
            matrix (np.array): Synteny matrix
        """
        total = len(self.clusters)
        matrix = np.zeros((total, total))
        for i, one in enumerate(self.clusters):
            for j, two in enumerate(self.clusters):
                if i == j:
                    continue
                matrix[i, j] = self.synteny(one, two, i=i)
        if normalise:
            # Explicitly check for zero so empty alignments don't cause issues
            max_value = matrix.max()
            if max_value > 0:
                matrix /= max_value
        if as_distance:
            maximum = 1 if normalise else matrix.max()
            matrix = maximum - matrix
            np.fill_diagonal(matrix, 0)
        return matrix

    def format_matrix(self, delimiter=",", **kwargs):
        """Generates a formatted distance matrix of stored clusters."""
        matrix = self.matrix(**kwargs)
        names = [self.clusters[uid].name for uid in self.clusters]
        header = ["", *names]
        result = [header]
        for index, row in enumerate(matrix):
            name = names[index]
            values = [str(value) for value in row]
            result.append([name, *values])
        return "\n".join(delimiter.join(row) for row in result)

    def order(self, i=0.5, method="ward"):
        """Determines optimal order of clusters using hierarchical clustering.

        When only a single cluster is stored, skips clustering and returns 0.
        """
        if len(self.clusters) == 1:
            return [0]
        if not self.alignments:
            return list(range(len(self.clusters)))
        matrix = self.matrix(i=i, normalise=True, as_distance=True)
        linkage = hierarchy.linkage(squareform(matrix), method=method)
        return hierarchy.leaves_list(linkage)[::-1]


class Alignment(Serializer):
    """An alignment between two gene clusters.

    Attributes:
        links (list): list of Gene-Gene 'links' (i.e. alignments)
    """

    def __init__(self, uid=None, query=None, target=None, links=None):
        self.uid = uid if uid else str(uuid.uuid4())
        self.query = query
        self.target = target
        self.links = links if links else []

    def to_dict(self, uids_only=False):
        return {
            "uid": self.uid,
            "query": self.query.uid if uids_only else self.query.to_dict(),
            "target": self.target.uid if uids_only else self.target.to_dict(),
            "links": [link.uid if uids_only else link.to_dict() for link in self.links]
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            query=load_child(d["query"], Cluster),
            target=load_child(d["target"], Cluster),
            links=load_children(d["links"], Link),
        )

    def __str__(self):
        return self.format()

    def format(
        self,
        decimals=4,
        delimiter=None,
        alignment_headers=True,
        link_headers=False,
    ):
        return format_alignment(
            self,
            decimals=decimals,
            delimiter=delimiter,
            alignment_headers=alignment_headers,
            link_headers=link_headers,
        )

    def contains(self, gene):
        """Return True if the given gene is in this cluster alignment."""
        return any(gene in (link.query, link.target) for link in self.links)

    @property
    def score(self):
        """Calculates the cumulative identity of this alignment."""
        if not self.links:
            raise ValueError("Alignment has no links")
        total = sum(link.identity for link in self.links)
        count = len(self.links)
        return total / count

    def add_link(self, query, target, identity, similarity):
        """Instantiate a new Link from a Gene alignment and save."""
        link = Link(
            query=query,
            target=target,
            identity=identity,
            similarity=similarity
        )
        self.links.append(link)


class Link(Serializer):
    """An alignment link between two Gene objects."""

    def __init__(self, uid=None, query=None, target=None, identity=None, similarity=None):
        self.uid = uid if uid else str(uuid.uuid4())
        self.query = query
        self.target = target
        self.identity = identity
        self.similarity = similarity

    def __str__(self):
        return self.format("\t")

    def values(self):
        return [self.query.name, self.target.name, self.identity, self.similarity]

    def to_dict(self, uids_only=False):
        return {
            "uid": self.uid,
            "query": self.query.uid if uids_only else self.query.to_dict(),
            "target": self.target.uid if uids_only else self.target.to_dict(),
            "identity": self.identity,
            "similarity": self.similarity,
        }

    @classmethod
    def from_dict(cls, d):
        d["query"] = load_child(d["query"], Gene)
        d["target"] = load_child(d["target"], Gene)
        return cls(**d)


class Group(Serializer):
    """A gene homology group."""

    def __init__(self, uid=None, label=None, genes=None, hidden=False, colour=None):
        self.uid = uid if uid else str(uuid.uuid4())
        self.label = label if label else f"Group {self.uid}"
        self.genes = genes if genes else []
        self.hidden = hidden
        self.colour = colour

    def to_dict(self):
        return {
            "uid": self.uid,
            "label": self.label,
            "genes": self.genes,
            "hidden": self.hidden,
            "colour": self.colour,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            uid=d.get("uid"),
            label=d.get("label"),
            genes=d.get("genes"),
            hidden=d.get("hidden"),
            colour=d.get("colour"),
        )
