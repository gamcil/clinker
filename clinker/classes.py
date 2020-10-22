#!/usr/bin/env python3

"""
Module to store classes used across the package.

Cameron Gilchrist
"""

import itertools
import warnings

from pathlib import Path

from Bio import SeqIO, SeqRecord, BiopythonParserWarning

# ignore malformed locus warnings
warnings.simplefilter('ignore', BiopythonParserWarning)


def find_qualifier(valid_values, qualifiers):
    """Try to find a valid value in SeqFeature.qualifiers"""
    for value in valid_values:
        if value in qualifiers:
            return qualifiers[value][0]
    return None


def parse_genbank(path):
    path = Path(path)
    with path.open() as fp:
        records = SeqIO.parse(fp, "genbank")
        cluster = Cluster.from_seqrecords(*records, name=path.stem)
    return cluster


def parse_files(paths):
    clusters = []
    for path in paths:
        if Path(path).is_dir():
            files = Path(path).glob("*")
            _clusters = parse_files(files)
            clusters.extend(_clusters)
        else:
            cluster = parse_genbank(path)
            clusters.append(cluster)
    return clusters


class Cluster:
    """The Cluster class stores Proteins

    Attributes:
        genes (list): list of Gene objects
        loci (list): index ranges of genomic loci corresponding
          to Gene objects in self.genes
    """

    id_iter = itertools.count()

    def __init__(self, name, loci, uid=None):
        self.uid = uid if uid else next(Cluster.id_iter)
        self.name = name
        self.loci = loci

    def to_dict(self):
        return {
            "uid": self.uid,
            "name": self.name,
            "loci": [locus.to_dict() for locus in self.loci]
        }

    @classmethod
    def from_seqrecords(cls, *args, name=None):
        """Instantiates a Cluster from BioPython SeqRecord/s.

        Arguments:
            args (SeqRecord): BioPython SeqRecord object/s
            mode (str): which sequence to extract ('protein' or 'gene')
        """
        loci = [Locus.from_seqrecord(record) for record in args]
        return cls(name=name if name else loci[0].name, loci=loci)

    def get_gene(self, name):
        for locus in self.loci:
            gene = locus.get_gene(name)
            if gene:
                return gene


class Locus:
    """A cluster locus."""

    id_iter = itertools.count()

    def __init__(self, name, genes, start=None, end=None, uid=None):
        self.uid = uid if uid else next(Locus.id_iter)
        self.name = name
        self.genes = genes
        self.start = start
        self.end = end

    def __str__(self):
        return f"{self.name}:{self.start}-{self.end}"

    def to_dict(self):
        return {
            "uid": self.uid,
            "name": self.name,
            "start": self.start,
            "end": self.end,
            "genes": [gene.to_dict() for gene in self.genes],
        }

    @classmethod
    def from_seqrecord(cls, record):
        """Builds a new Locus from a BioPython SeqRecord."""
        if not isinstance(record, SeqRecord.SeqRecord):
            raise NotImplementedError('Supplied argument is not a valid SeqRecord object')
        genes = [
            Gene.from_seqfeature(feature, record)
            for feature in record.features
            if feature.type == "CDS"
        ]
        return cls(name=record.name, start=0, end=len(record), genes=genes)

    def get_gene(self, name):
        for gene in self.genes:
            if gene.name == name:
                return gene


class Gene:
    """Location, annotation and attachment points for drawing links."""

    id_iter = itertools.count()

    def __init__(
        self,
        name=None,
        start=None,
        end=None,
        strand=None,
        function=None,
        sequence=None,
        translation=None,
        uid=None,
    ):
        self.uid = uid if uid else next(Gene.id_iter)
        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence
        self.translation = translation
        self.smcog = {}
        self.function = function if function else 'Hypothetical'

    def to_dict(self):
        return {
            "uid": self.uid,
            "name": self.name,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "sequence": self.sequence,
            "translation": self.translation,
            "smcog": self.smcog,
            "function": self.function,
        }

    @property
    def ordered_smcogs(self):
        """Returns an ordered list of smCOG hits."""
        ordered = sorted(self.smcog, key=lambda x: self.smcog[x][0])
        return [(key, *self.smcog[key]) for key in ordered]

    @classmethod
    def from_seqfeature(cls, feature, record):
        """Builds a new Gene object from a BioPython SeqFeature.

        Parameters:
            feature (SeqFeature): BioPython SeqFeature object
            record (SeqRecord): BioPython SeqRecord object (parent of feature)
        """
        name = find_qualifier(['protein_id', 'locus_tag', 'ID'], feature.qualifiers)
        if not name:
            raise ValueError(
                "Could not determine a valid identifier"
                f" from a CDS SeqFeature in {record.id}"
            )
        sequence = feature.extract(record.seq)
        translation = find_qualifier(
            ["translation"],
            feature.qualifiers,
        )
        if not translation:
            translation = sequence.translate()
        return cls(
            name=name,
            sequence=str(sequence),
            translation=str(translation),
            start=int(feature.location.start),
            end=int(feature.location.end),
            strand=feature.location.strand,
        )
