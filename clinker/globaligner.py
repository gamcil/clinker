#!/usr/bin/env python3

"""
Align protein sequences in Clusters or lists of Proteins.

Cameron Gilchrist
"""

from collections import defaultdict, OrderedDict, namedtuple
from itertools import combinations, permutations

from Bio import Align
from Bio.Align import substitution_matrices

from classes import Cluster


class ClusterAligner:
    """ Performs and stores alignments.

        Attributes:
             _clusters (dict): internal ordered dict of Cluster objects stored
                               by their name attribute
           _alignments (list): internal list of Alignment objects
   _aligner (PairwiseAligner): BioPython Align.PairwiseAligner object used
                               for aligning sequences
    _alignment_indices (dict): indices of Alignment objects in _alignments
                               stored using Cluster.name attributes as keys
        _cluster_names (dict): tuples of Cluster.name attributes stored using
                               _alignment indices as keys
    """

    aligner_default = {
        "mode": "global",
        "substitution_matrix": substitution_matrices.load("BLOSUM62"),
        "open_gap_score": -10,
        "extend_gap_score": -0.5,
    }

    def __init__(self, aligner_config=None):
        self._alignments = []
        self._clusters = OrderedDict()
        self._alignment_indices = defaultdict(dict)
        self._cluster_names = defaultdict(dict)
        self._aligner = Align.PairwiseAligner()

        if aligner_config:
            self.configure_aligner(**aligner_config)
        else:
            self.configure_aligner(**self.aligner_default)

    def __repr__(self):
        """Print all alignments currently stored in the instance."""
        if not self._alignments:
            return f"No alignments are currently stored."

        strings = set()
        for parent, children in self._alignment_indices.items():
            for child in children:
                alignment = str(self.get_alignment(parent, child))

                if alignment not in strings:
                    header = f"{parent} vs {child}"
                    separator = "-" * len(header)

                    strings.add(f"{header}\n{separator}\n{alignment}")

        return "\n\n".join(strings)

    def add_clusters(self, *clusters):
        """ Add new Cluster object/s to the ClusterAligner.

            Arguments:
                *clusters (list): variable number of Cluster objects
        """
        for cluster in clusters:
            if not isinstance(cluster, Cluster):
                raise ValueError("Supplied cluster is not a valid Cluster object")
            if cluster.name not in self._clusters:
                self._clusters[cluster.name] = cluster

    def align_stored_clusters(self, cutoff=0.3):
        """ Run alignments on clusters stored in the ClusterAligner.
        """
        for one, two in combinations(self._clusters.values(), 2):

            if (
                one.name in self._alignment_indices[two.name]
                or two.name in self._alignment_indices[one.name]
            ):
                # alignment already performed, skip
                continue

            self.add_alignment(
                one, two, align_two_clusters(one, two, self._aligner, cutoff=cutoff)
            )

    def configure_aligner(self, **kwargs):
        """Change properties on the BioPython.PairwiseAligner object.

        This function takes any keyword argument and assumes
        they correspond to valid properties on the PairwiseAligner.
        Refer to BioPython documentation for these.
        """
        valid_attributes = set(dir(self._aligner))

        for key, value in kwargs.items():
            if key not in valid_attributes:
                raise ValueError(
                    f'"{key}" is not a valid attribute of the BioPython'
                    "Align.PairwiseAligner class"
                )
            setattr(self._aligner, key, value)

    @property
    def aligner_settings(self):
        """ Return str representation of PairwiseAligner object, which
            lists all current settings.
        """
        return str(self._aligner)

    def form_alignment_string(self, index):
        """ Return a string representation of a stored Alignment.
        """
        one, two = self._cluster_names[index]

        header = f"{one} vs {two}"
        separator = "-" * len(header)
        alignment = self._alignments[index]

        return f"{header}\n{separator}\n{alignment}"

    def add_alignment(self, query, target, alignment):
        """ Add a new cluster alignment to the ClusterAligner.

            self._alignment_indices allows for Alignment indices to be
            retrieved from cluster names, regardless of order.

            self._cluster_names allows for Cluster names to be retrieved
            given the index of an Alignment in self._alignments
        """
        # Save Cluster object if not already stored
        self.add_clusters(query, target)

        # Update mapping dictionaries and save Alignment
        index = len(self._alignments)
        self._alignment_indices[query.name][target.name] = index
        self._alignment_indices[target.name][query.name] = index
        self._cluster_names[index] = (query.name, target.name)
        self._alignments.append(alignment)

    def get_alignment(self, one, two):
        """ Retrieve a cluster alignment.

            Arguments:
                one (str): name corresponding to a Cluster object
                two (str): name corresponding to a Cluster object
            Returns:
                Alignment object for the specified Clusters
        """
        try:
            # Since alignment indices are saved in both directions,
            # only need to check one way
            return self._alignments[self._alignment_indices[one][two]]
        except KeyError as exc:
            raise KeyError(
                "No alignment could be retrieved for these clusters"
            ) from exc

    def find_largest_cluster(self):
        """ Find the largest cluster stored in the ClusterAligner.

            This is used for scaling purposes - i.e. if width=1000, pixel 0
            will be the start of the first loci of this cluster, and pixel 1000
            will be the end of the last loci.
        """
        largest_length, largest_cluster = 0, None
        for cluster in self._clusters.values():
            if cluster.length > largest_length:
                largest_length = cluster.length
                largest_cluster = cluster
        return largest_cluster

    def get_cluster(self, cluster):
        """ Retrieve a Cluster object given its name.

            Arguments:
                cluster (str): name attribute of a Cluster object
        """
        try:
            return self._clusters[cluster]
        except KeyError as exc:
            raise KeyError("No cluster with this name could be retrieved") from exc

    # TODO: check for number of genes vs pure identity score
    def compute_display_order(self, order="first_as_seed"):
        """ Return sorted cluster name list based on specified order scheme.

            Arguments:
                order (int):
                    'maximize'
                        Return the permutation of clusters that maximizes
                        the sum of alignment scores between each cluster pair.

                    'first_as_seed'
                        Use first added cluster as seed (i.e. first element),
                        and maximize scores for subsequent clusters.

                    'added'
                        Returns clusters in the order in which they were added.

            Returns:
                permutation (list): ordered list of Cluster names
        """
        total_clusters = len(self._clusters)

        if order == "maximize":
            # Brute force; scores of all possible permutations
            best_score = 0
            best_permutation = None
            for permutation in permutations(self._clusters.keys()):
                score = 0
                for i in range(1, total_clusters):
                    try:
                        score += self.get_alignment(
                            permutation[i - 1], permutation[i]
                        ).score
                    except ValueError:
                        # no links
                        continue
                if score > best_score:
                    best_score = score
                    best_permutation = permutation

            return list(best_permutation)

        if order == "first_as_seed":
            # First cluster as seed, then best remaining cluster at each
            # subsequent step
            seed, *clusters = list(self._clusters.keys())
            permutation = [seed]

            for i in range(total_clusters - 1):
                max_score = 0
                max_cluster = None
                for target in clusters:
                    # Skip if already in the permutation
                    if target in permutation:
                        continue

                    # Grab alignment between last added cluster and next
                    try:
                        alignment = self.get_alignment(permutation[-1], target)
                    except KeyError:
                        continue
                    try:
                        if alignment.score > max_score:
                            max_score = alignment.score
                            max_cluster = target
                    except ValueError:
                        continue

                if max_cluster:
                    permutation.append(max_cluster)

            return permutation

        if order == "added":
            # Since an OrderedDict is used internally, just return keys
            return list(self._clusters.keys())

        raise ValueError(
            "Value for order argument must be either"
            ' "maximize", "first_as_seed", or "added"'
        )


class Alignment:
    """An alignment between two gene clusters.

    An instance of this class is agnostic to the clusters
    that have been aligned to create it; it will typically
    be stored and accessed via the ClusterAligner class.

    Attributes:
        links (list): list of Gene-Gene 'links' (i.e. alignments)
    """

    # Use a namedtuple to store gene-gene homologies
    Link = namedtuple("Link", ["query", "target", "identity", "similarity"])

    def __init__(self, links=None):
        self.links = links if links else []

    def __str__(self):
        return "\n".join(
            f"{link.query},"
            f"{link.target},"
            f"{link.identity:.4f},"
            f"{link.similarity:.4f}"
            for link in self.links
        )

    def find_gene(self, gene_name):
        """ Return True if the given gene is in this cluster alignment.
        """
        for link in self.links:
            if gene_name in (link.query, link.target):
                return True
        return False

    @property
    def score(self):
        """Calculate 'score' of this alignment by summing identities."""
        if not self.links:
            raise ValueError("Alignment has no links")
        return sum(link.identity + link.similarity for link in self.links) / len(
            self.links
        )

    def add_link(self, query, target, identity, similarity):
        """Instantiate a new Link from a Gene alignment and save."""
        self.links.append(self.Link(query, target, identity, similarity))


def align_clusters(*args, cutoff=0.3, aligner_config=None):
    """Convenience function for directly aligning Cluster object/s.

    Initialises a ClusterAligner, adds Cluster/s, then runs alignments
    and returns the ClusterAligner.

    Args:
        *args: Cluster or list of Protein objects
        aligner_config (dict): keyword arguments to use when setting
                               up the BioPython.PairwiseAligner object
        cutoff (float): decimal identity cutoff for saving an alignment
    Returns:
        aligner (ClusterAligner): instance of ClusterAligner class which
                                  contains all cluster alignments
    e.g.
        align_sequence_groups(cluster1, cluster2, ..., clusterN)
    """
    if len(args) < 2:
        raise ValueError("Must provide 2 or more clusters")
    aligner = ClusterAligner(aligner_config)
    aligner.add_clusters(*args)
    aligner.align_stored_clusters(cutoff)
    return aligner


def align_two_clusters(one, two, aligner, cutoff=0.3):
    """Perform pairwise global alignments between two Clusters.

    Args:
        one (Cluster): Cluster instance
        two (Cluster): Cluster instance
    aligner (Align.PairwiseAligner): BioPython aligner object
         cutoff (int): identity threshold for saving a link

    Returns:
        homologies (list): identity and similarity decimals
                           for each pairwise comparison that
                           met the identity cutoff
        score (int): sum of identities in homologies, for use
                     as a homology score
    """

    def compute_identity(alignment):
        """ Calculate sequence identity and similarity of
            an alignment.
        """
        # Aligned strings aren't stored separately, have to split
        one, _, two, _ = str(alignment).split("\n")
        length = len(one)

        # Amino acid similarity groups
        similar_acids = [
            {"G", "A", "V", "L", "I"},
            {"F", "Y", "W"},
            {"C", "M"},
            {"S", "T"},
            {"K", "R", "H"},
            {"D", "E", "N", "Q"},
            {"P"},
        ]

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

    alignment = Alignment()
    for gene_one in one.genes:
        for gene_two in two.genes:

            # Run the alignment, compute identity and similarity
            identity, similarity = compute_identity(
                aligner.align(gene_one.sequence, gene_two.sequence)[0]
            )

            # Save if it meets the threshold
            if identity >= cutoff:
                alignment.add_link(gene_one.name, gene_two.name, identity, similarity)
    return alignment
