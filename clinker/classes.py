#!/usr/bin/env python3

"""
Module to store classes used across the package.

Cameron Gilchrist
"""

import json
import warnings
import uuid

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


def subdict(d, keys):
    """Creates a sub-dictionary of a parent dictionary given set of keys."""
    funcs = ("lower", "upper", "title")
    sub = {}
    for key, value in d.items():
        if key in keys or any(getattr(key, func) in keys for func in funcs):
            if isinstance(value, list):
                sub[key] = value[0]
            else:
                sub[key] = value
    return sub


def get_value(d, keys):
    for key in keys:
        if key in d:
            return d[key]


def parse_genbank(path):
    path = Path(path)
    with path.open() as fp:
        records = SeqIO.parse(fp, "genbank")
        cluster = Cluster.from_seqrecords(*records, name=path.stem)
    return cluster


def find_files(paths, recurse=True, level=0):
    files = []
    for path in paths:
        if Path(path).is_dir():
            if level == 0 or recurse:
                new = Path(path).glob("*")
                _files = find_files(new, recurse=recurse, level=level + 1)
                files.extend(_files)
        else:
            files.append(path)
    return files


def parse_files(paths):
    return [parse_genbank(path) for path in paths]


def get_children(children, uids_only=False):
    return [
        child.uid
        if uids_only
        else child.to_dict()
        for child in children
    ]


def load_child(child, thing):
    if not hasattr(thing, "from_dict"):
        raise NotImplementedError
    return child if isinstance(child, str) else thing.from_dict(child)


def load_children(children, thing):
    return [load_child(child, thing) for child in children]


class Serializer:
    """JSON serialisation mixin class.

    Classes that inherit from this class should implement `to_dict` and
    `from_dict` methods.
    """

    def to_dict(self):
        """Serialises class to dict."""
        raise NotImplementedError

    @classmethod
    def from_dict(self, d):
        """Loads class from dict."""
        raise NotImplementedError

    def to_json(self, fp=None, **kwargs):
        """Serialises class to JSON."""
        d = self.to_dict()
        if fp:
            json.dump(d, fp, **kwargs)
        else:
            return json.dumps(d, **kwargs)

    @classmethod
    def from_json(cls, js):
        """Instantiates class from JSON handle."""
        if isinstance(js, str):
            d = json.loads(js)
        else:
            d = json.load(js)
        return cls.from_dict(d)


class Cluster(Serializer):
    """The Cluster class stores Proteins

    Attributes:
        genes (list): list of Gene objects
        loci (list): index ranges of genomic loci corresponding
          to Gene objects in self.genes
    """

    def __init__(self, name, loci, uid=None):
        self.uid = uid if uid else str(uuid.uuid4())
        self.name = name
        self.loci = loci

    def to_dict(self, uids_only=False):
        return {
            "uid": self.uid,
            "name": self.name,
            "loci": get_children(self.loci, uids_only=uids_only)
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["name"],
            load_children(d["loci"], Locus),
            uid=d.get("uid")
        )

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


class Locus(Serializer):
    """A cluster locus."""

    def __init__(self, name, genes, start=None, end=None, uid=None):
        self.uid = uid if uid else str(uuid.uuid4())
        self.name = name
        self.genes = genes
        self.start = start
        self.end = end

    def __str__(self):
        return f"{self.name}:{self.start}-{self.end}"

    def to_dict(self, uids_only=False):
        return {
            "uid": self.uid,
            "name": self.name,
            "start": self.start,
            "end": self.end,
            "genes": get_children(self.genes, uids_only=uids_only),
        }

    @classmethod
    def from_dict(cls, d):
        return cls(
            d["name"],
            load_children(d["genes"], Gene),
            start=d["start"],
            end=d["end"],
            uid=d.get("uid"),
        )

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
        genes = [gene for gene in genes if gene]
        return cls(name=record.name, start=0, end=len(record), genes=genes)

    def get_gene(self, name):
        for gene in self.genes:
            if gene.name == name:
                return gene


class Gene(Serializer):
    """Location, annotation and attachment points for drawing links."""

    def __init__(
        self,
        uid=None,
        label=None,
        names=None,
        start=None,
        end=None,
        strand=None,
        sequence=None,
        translation=None,
    ):
        self.uid = uid if uid else str(uuid.uuid4())
        self.label = label if label else self.uid
        self.names = names if names else {}
        self.start = start
        self.end = end
        self.strand = strand
        self.sequence = sequence
        self.translation = translation

    def to_dict(self):
        return {
            "uid": self.uid,
            "label": self.label,
            "names": self.names,
            "start": self.start,
            "end": self.end,
            "strand": self.strand,
            "sequence": self.sequence,
            "translation": self.translation,
        }

    @classmethod
    def from_dict(cls, d):
        return cls(**d)

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
        if "pseudo" in feature.qualifiers:
            return
        tags = ("protein_id", "locus_tag", "id", "gene", "label", "name")
        names = subdict(feature.qualifiers, tags)
        sequence = feature.extract(record.seq)
        translation = find_qualifier(
            ["translation"],
            feature.qualifiers,
        )
        if not translation:
            translation = sequence.translate()
        return cls(
            names=names,
            label=get_value(names, tags),
            sequence=str(sequence),
            translation=str(translation),
            start=int(feature.location.start),
            end=int(feature.location.end),
            strand=feature.location.strand,
        )
