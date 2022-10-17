#!/usr/bin/env python3

"""
Module to store classes used across the package.

Cameron Gilchrist
"""

import json
import logging
import warnings
import uuid
import copy
import operator

from pathlib import Path

import gffutils
from gffutils import biopython_integration

from Bio import SeqIO, SeqRecord, BiopythonParserWarning
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation

# ignore malformed locus warnings
warnings.simplefilter('ignore', BiopythonParserWarning)

LOG = logging.getLogger(__name__)

FASTA_SUFFIXES = (".fa", ".fsa", ".fna", ".fasta", ".faa")
GBK_SUFFIXES = (".gbk", ".gb", ".genbank", ".gbf", ".gbff")
GFF_SUFFIXES = (".gtf", ".gff", ".gff3")


def merge_locations(parts):
    new_parts = []
    for i in range(0, len(parts), 2):
        a, b = parts[i], parts[i + 1]
        if a.end == b.start:
            c = FeatureLocation(a.start, b.end, a.strand, a.ref, a.ref_db)
            new_parts.append(c)
        else:
            new_parts.extend(new_parts[i: i + 1])
    return new_parts


def shift_origin(record: SeqRecord) -> SeqRecord:
    """Sets a new origin on a circular SeqRecord and shifts relevant SeqFeatures.

    Adapted from PyDNA
    https://github.com/BjornFJohansson/pydna/blob/75fa7e2c783eabfa32f580d7c0ea2620f907390b/src/pydna/dseqrecord.py#L1019
    """

    # Collect features going across the current origin
    spanning_origin = []
    for feature in record.features:
        if feature.type.lower() not in ('gene', 'cds'):
            continue
        if feature.location.parts[0].start < feature.location.parts[-1].end:
            continue
        spanning_origin.append(feature)

    if not spanning_origin:
        LOG.info("No features wrapping the origin")
        return record

    # New origin will just be before the start of the first feature spanning the origin
    shift = spanning_origin[0].location.parts[0].start

    LOG.info(f"Found {len(spanning_origin)} features wrapping origin")
    LOG.info(f"Setting new origin to {shift} (start of {feature.qualifiers.get('locus_tag')[0]})")

    # BioPython-recommended way of shifting origin, but deletes SeqFeatures which span the origin
    record = record[shift:] + record[:shift]
    length = len(record)

    if not shift % length:
        return record  # shift is a multiple of ln or 0

    shift = length - (shift % length)

    for feature in spanning_origin:
        new_parts = []
        for part in feature.location.parts:
            part += shift
            new_start = part.start % length
            new_end = part.end % length

            if new_start < new_end:
                fl = FeatureLocation(new_start, new_end, part.strand, part.ref, part.ref_db)
                new_parts.append(fl)

            # Fix other features if they now wrap the origin
            elif new_start > new_end:
                if part.strand == 1:
                    fls = [
                        FeatureLocation(new_start, length, part.strand, part.ref, part.ref_db),
                        FeatureLocation(0, new_end, part.strand, part.ref, part.ref_db),
                    ]
                else:
                    fls = [
                        FeatureLocation(0, new_end, part.strand, part.ref, part.ref_db),
                        FeatureLocation(new_start, length, part.strand, part.ref, part.ref_db),
                    ]
                new_parts.extend(fls)

        # Merge parts with identical end/start positions
        if len(new_parts) > 1:
            new_parts = merge_locations(new_parts)
        if len(new_parts) > 1:
            feature.location = CompoundLocation(final)
        else:
            feature.location = new_parts[0]

    record.features.extend(spanning_origin)
    return record


def find_fasta(gff_path):
    path = None
    for suffix in FASTA_SUFFIXES:
        _path = Path(gff_path).with_suffix(suffix)
        if _path.exists():
            path = _path
            break
    return path


def parse_fasta(path):
    with open(path) as fp:
        records = list(SeqIO.parse(fp, "fasta"))
    return records


def find_regions(directives):
    """Looks for ##sequence-region directives in a list of GFF3 directives."""
    regions = {}
    for directive in directives:
        if directive.startswith("sequence-region"):
            _, accession, start, end = directive.split(" ")
            regions[accession] = (int(start), int(end))
    return regions


def cluster_from_gff(path, ranges=None):
    """Parses a GFF3 file using GFFUtils."""

    # Check for FASTA file
    fasta_path = find_fasta(path)
    if not fasta_path:
        raise FileNotFoundError("Could not find matching FASTA file")

    # Parse FASTA and create GFFUtils database
    fasta = parse_fasta(fasta_path)
    gff = gffutils.create_db(
        str(path),
        ":memory:",
        force=True,
        merge_strategy="create_unique",
        sort_attribute_values=True
    )
    regions = find_regions(gff.directives)

    # Find features for each record in the FASTA file
    loci = []
    for record in fasta:
        # Check for matching ##sequence-region directive
        try:
            record_start, record_end = regions[record.id]
        except KeyError:
            record_start, record_end = 1, len(record)

        # Check for user-specified range
        if ranges and record.id in ranges:
            record_start, record_end = ranges[record.id]
            LOG.info("    Parsing range %s:%i-%i", record.id, record_start, record_end)

            # Adjust FASTA record to match record_start and record_end
            # -- Default: 0 to end of record
            # -- ##sequence-region: start to end of directive
            # -- User-specified range: start to end of range
            record = record[record_start - 1: record_end]

        # Zero-index the start of the record
        record_start -= 1

        # Extract features from record within range
        region = gff.region(
            seqid=record.id,
            featuretype="CDS",
            start=record_start,
            end=record_end,
            completely_within=True,
        )
        features = sorted(region, key=lambda f: f.start)
        if not features:
            LOG.warning("Found no CDS features in %s [%s]", record.id, path)

        previous = None
        for feature in features:
            # Check if this feature is part of the previous one for merging
            seqid = feature.attributes["ID"][0]
            same_feature = previous == seqid
            if not previous:
                previous = seqid

            # Normalise Feature location based on ##sequence-region directive.
            # Necessary for extracted GFF3 files that still store coordinates
            # relative to the entire region. If no sequence-region directive
            # is found, assumes 1 (i.e. default sequence start).
            # Note: to_seqfeature automatically zero indexes coordinates, which
            # does not happen by default in GFFUtils, hence no -1 here
            feature = biopython_integration.to_seqfeature(feature)
            feature.location = FeatureLocation(
                feature.location.start - record_start,
                feature.location.end - record_start,
                strand=feature.location.strand
            )

            # Either merge with previous feature, or append it
            if same_feature:
                if feature.location.strand == 1:
                    record.features[-1].location += feature.location
                else:
                    # Must be in biological order
                    old, new = record.features[-1].location, feature.location
                    record.features[-1].location = new + old
            else:
                record.features.append(feature)
                previous = seqid

        # Try to trace back from CDS to parent gene feature for actual
        # gene coordinates. If not found (e.g. malformed GFF without ID= and parent=
        # features), warns user and defaults to CDS start/end.
        genes = []
        for feature in record.features:
            parents = [p for p in gff.parents(gff[feature.id], featuretype="gene")]
            if parents:
                parent, *_ = parents
                # e.g. CDS is within range, but gene UTR is not
                if parent.start < record_start or parent.end > record_end:
                    continue
                start = parent.start
                end = parent.end
            else:
                LOG.warning(
                    f"Could not find parent gene of {feature.id}."
                    " Using coding sequence coordinates instead."
                )
                start = feature.location.start + record_start
                end = feature.location.end + record_start
            gene = Gene.from_seqfeature(feature, record, start=start, end=end)
            genes.append(gene)

        locus = Locus(record.id, genes, start=record_start, end=record_end)
        loci.append(locus)

    return Cluster(Path(path).stem, loci)


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


def cluster_from_genbank(path, ranges=None, set_origin=True):
    """Instantiates a new Cluster object from a GenBank file at `path`.

    Args:
        path (str): Path to GenBank file
        ranges (dict): Dictionary of scaffold ranges
    Returns:
        Cluster
    """
    path = Path(path)
    loci = []
    with path.open() as fp:
        for record in SeqIO.parse(fp, "genbank"):
            if ranges and record.id in ranges:
                start, end = ranges[record.id]
                locus = Locus.from_seqrecord(record, start=start - 1, end=end, set_origin=set_origin)
            else:
                locus = Locus.from_seqrecord(record, set_origin=set_origin)
            loci.append(locus)
    return Cluster(path.stem, loci)


def find_files(paths, recurse=True, level=0):
    files = []
    for path in paths:
        _path = Path(path)
        if _path.is_dir():
            if level == 0 or recurse:
                new = Path(path).glob("*")
                _files = find_files(new, recurse=recurse, level=level + 1)
                files.extend(_files)
        else:
            if _path.exists() and _path.suffix.lower() in GBK_SUFFIXES + GFF_SUFFIXES:
                files.append(path)
    return files


def parse_files(paths, ranges=None, set_origin=True):
    clusters = []
    for path in paths:
        LOG.info("  %s", Path(path).name)
        if Path(path).suffix.lower() in GBK_SUFFIXES:
            cluster = cluster_from_genbank(path, ranges=ranges, set_origin=set_origin)
        elif Path(path).suffix.lower() in GFF_SUFFIXES:
            cluster = cluster_from_gff(path, ranges=ranges)
        else:
            raise TypeError("File %s does not have GenBank or GFF3 extension")
        clusters.append(cluster)
    return clusters


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


def find_overlapping_location(feature, locations):
    for location in locations:
        if (
            feature.location.start >= location.start
            and feature.location.end <= location.end
        ):
            return location


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
        self.loci = loci if loci else []

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
        """
        loci = [Locus.from_seqrecord(record) for record in args]
        return cls(name=name if name else loci[0].name, loci=loci)

    def get_gene(self, label):
        for locus in self.loci:
            gene = locus.get_gene(label)
            if gene:
                return gene


class Locus(Serializer):
    """A cluster locus."""

    def __init__(
        self,
        name,
        genes,
        start=None,
        end=None,
        uid=None,
    ):
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
    def from_seqrecord(cls, record, set_origin=True, start=None, end=None):
        """Builds a new Locus from a BioPython SeqRecord."""

        if not isinstance(record, SeqRecord.SeqRecord):
            raise TypeError("Expected SeqRecord object")

        if set_origin and record.annotations.get('topology') == 'circular':
            LOG.warning(f"{record.id} is circular, checking for genes spanning the origin")
            LOG.warning("To disable this behaviour, use --dont_set_origin")
            record = shift_origin(record)

        ranged = start or end
        start = start if start else 0
        end = end if end else len(record)
        
        # Manually filter SeqRecord if start or end is specified.
        # Can directly slice SeqRecord objects, but will automatically adjust
        # feature locations - we want them as is, since they'll be normalised
        # by the range during visualisation.
        if ranged:
            LOG.info("    Parsing range %s:%i-%i", record.id, start, end)
            record.features = [
                feature
                for feature in record.features
                if feature.location.start >= start and feature.location.end <= end
            ]
            record.seq = record.seq[start - 1: end]

        # Find all CDS SeqFeature and gene FeatureLocations
        features = [f for f in record.features if f.type == "CDS"]
        locations = [f.location for f in record.features if f.type == "gene"]

        # Trace CDS features back to genes for real locations, create Gene objects
        genes = []
        for feature in features:
            match = find_overlapping_location(feature, locations)
            if not match:
                LOG.warning(
                    f"Could not find parent gene of {feature.id}."
                    " Using coding sequence coordinates instead."
                )
            gene = Gene.from_seqfeature(
                feature,
                record,
                start=int(match.start) if match else None,
                end=int(match.end) if match else None,
            )
            if gene:
                genes.append(gene)

        # Instantiate the Locus object
        return cls(
            name=record.name,
            start=start if start else 0,
            end=end if end else len(record),
            genes=genes,
        )

    def get_gene(self, label):
        for gene in self.genes:
            if gene.label == label:
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
    def from_seqfeature(
        cls,
        feature,
        record,
        start=None,
        end=None,
    ):
        """Builds a new Gene object from a BioPython SeqFeature.

        Parameters:
            feature (SeqFeature): BioPython SeqFeature object
            record (SeqRecord): BioPython SeqRecord object (parent of feature)
        """
        tags = ("protein_id", "locus_tag", "id", "gene", "label", "name")
        qualifiers = {
            k: v[0] if isinstance(v, list) else v
            for k, v in feature.qualifiers.items()
        }
        sequence = feature.extract(record.seq)
        translation = qualifiers.pop("translation", sequence.translate())
        return cls(
            names=qualifiers,
            label=get_value(qualifiers, tags),
            sequence=str(sequence),
            translation=str(translation),
            start=start if start else int(feature.location.start),
            end=end if end else int(feature.location.end),
            strand=feature.location.strand,
        )
