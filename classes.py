#!/usr/bin/env python3

"""
Module to store classes used across the package.

Cameron Gilchrist
"""

import warnings

from collections import OrderedDict

from Bio import SeqRecord, BiopythonParserWarning

# ignore malformed locus warnings
warnings.simplefilter('ignore', BiopythonParserWarning)


class Cluster:
    """ The Cluster class stores Proteins

        Attributes:
            genes (list): list of Gene objects
            loci (list): index ranges of genomic loci corresponding
                         to Gene objects in self.genes
    """
    def __init__(self, name, genes, loci):
        self.name = name
        self.loci = loci
        self.genes = []
        self.length = 0

        # Sort genes by start location, but ensure separate loci
        # remain separate in this list
        for start, end in self.loci.values():
            self.genes.extend(
                sorted(genes[start:end],
                       key=lambda gene: gene.location[0])
            )

            # Keep track of total cluster length, minus locus spacing
            self.length += (self.genes[end - 1].location[1] -
                            self.genes[start].location[0])

        # Start and end Gene indices for drawing purposes
        # corresponding to homology regions
        self.drawing_bounds = [None, None]

    def spaced_length(self, locus_spacing=40):
        """ Return the total length of this cluster, including locus spacing
        """
        return self.length + locus_spacing * len(self.loci)

    def map_cluster(self, locus_spacing=40, scale_factor=1.0):
        """ Generate a scaled map of the genes in the loci of this Cluster.

            Each gene is scaled to the locus it is on, and saved
            as a percentage of the entire locus for scaling purposes.

            e.g.
            'Locus1' : {
                'length': 15000,
                'genes': {
                    'GENE_0001': (10, 100, -1),
                    'GENE_0002': (200, 300, 1),
                    ...
                }
            },
            'Locus2': {
                ...
            }

            This can then be passed to the visualizer module.
        """
        loci = {}
        prev = None
        for locus, (start, end) in self.loci.items():
            # Calculate width of this loci
            locus_length = (self.genes[end - 1].location[1] -
                            self.genes[start].location[0])

            # If multiple loci, add buffer between them
            horizontal_offset = (loci[prev]['length'] + locus_spacing
                                 if prev else 0)
            loci[locus] = {
                'start': horizontal_offset,
                'end': horizontal_offset + locus_length*scale_factor,
                'genes': OrderedDict()
            }

            # Save scaled start/end; want to be able to separately control
            # vertical spacing and gene arrow shape (arrow tip size etc)
            for gene in self.genes[start:end]:
                loci[locus]['genes'][gene.name] = (
                    horizontal_offset + gene.location[0]*scale_factor,
                    horizontal_offset + gene.location[1]*scale_factor,
                    gene.location[2],
                    gene.function
                )
            prev = locus

        return loci

    @classmethod
    def from_seqrecords(cls, *args, mode='protein'):
        """ Instantiate a Cluster from BioPython SeqRecord/s.

            Handles multiple records as separate genomic loci.
            If multiple are provided, this cluster will be drawn
            separated on the same line.

            Arguments:
               *args (SeqRecord): BioPython SeqRecord object/s
                mode (str): type of sequence to extract and
                            compare, either 'protein' or 'gene'
        """
        genes = []
        loci = {}
        for record in args:
            if not isinstance(record, SeqRecord.SeqRecord):
                raise ValueError(
                    'Supplied argument is not a valid SeqRecord object'
                )

            name = record.name

            feature_count = 0
            for feature in record.features:
                if mode == 'gene' and feature.type == 'gene':
                    gene = Gene.from_gene_seqfeature(feature, record)

                elif mode == 'protein' and feature.type == 'CDS':
                    gene = Gene.from_cds_seqfeature(feature, record)

                else:
                    continue

                genes.append(gene)
                feature_count += 1

            # Save the indices of the genes at this locus
            start = len(loci)
            end = start + feature_count
            loci[record.id] = (start, end)

        return cls(name=name, genes=genes, loci=loci)


def find_qualifier(valid_values, qualifiers):
    """ Try to find a valid value in SeqFeature.qualifiers
    """
    for value in valid_values:
        if value in qualifiers:
            return qualifiers[value][0]
    return None


class Gene:
    """ Location, annotation and attachment points for drawing links.
    """
    def __init__(self, name=None, location=None,
                 function=None, sequence=None):
        self.name = name
        self.location = location
        self.function = function if function else 'Hypothetical'
        self.sequence = sequence
        self.attachments = [None, None]  # start, end
        self.smcog = {}

    @property
    def ordered_smcogs(self):
        """ Return an ordered list of smCOG hits.
        """
        ordered = sorted(self.smcog, key=lambda x: self.smcog[x][0])
        return [(key, *self.smcog[key]) for key in ordered]

    @classmethod
    def from_gene_seqfeature(cls, feature, record):
        """ Instantiate a Gene object from a gene SeqFeature.
        """
        # Try and get an identifier
        name = find_qualifier(['locus_tag', 'ID'],
                              feature.qualifiers)
        if not name:
            raise ValueError(
                'Could not determine a valid identifier'
                f' from a gene SeqFeature in {record.id}'
            )

        return cls(
            name=name,
            location=(feature.location.start,
                      feature.location.end,
                      feature.location.strand),
            sequence=feature.extract(record.seq)
        )

    @classmethod
    def from_cds_seqfeature(cls, feature, record):
        """ Instantiate a Gene object from a CDS SeqFeature.
        """
        name = find_qualifier(['protein_id', 'locus_tag', 'ID'],
                              feature.qualifiers)
        if not name:
            raise ValueError(
                'Could not determine a valid identifier'
                f' from a CDS SeqFeature in {record.id}'
            )

        translation = find_qualifier(['translation'], feature.qualifiers)
        if not translation:  # manually translate with BioPython
            translation = feature.extract(record.seq).translate()

        return cls(
            name=name,
            location=(feature.location.start,
                      feature.location.end,
                      feature.location.strand),
            sequence=translation
        )
