#!/usr/bin/env python3

"""
Module for annotation of Gene objects.

Cameron Gilchrist
"""

import shutil
import subprocess

from collections import defaultdict
from tempfile import NamedTemporaryFile

HMMSEARCH_PATH = shutil.which('hmmsearch')
HMMSEARCH_DATA = 'smcogs/smcog.hmm'

if not HMMSEARCH_PATH:
    raise FileNotFoundError(
        'Could not find hmmsearch on your system PATH!'
        ' Please ensure it is installed.'
    )


def assign_function(domains):
    """ Assign a function to a Gene object given its annotations.

        Arguments:
            domains (list): all domain hits for a Gene
        Returns:
            function (str): assigned biosynthetic function based on
                            smCOG profile hits
    """
    if not domains:
        return 'Hypothetical'

    # Polyketide synthase (PKS), just check PKS+AT
    pks = {'SMCOG1022', 'SMCOG1021'}

    # Non-ribosomal peptide synthase (NRPS), just check A
    nrps = {'SMCOG1002'}

    # PKS-NRPS hybrids, check for PKS+AT+A
    hybrid = {'SMCOG1022', 'SMCOG1141', 'SMCOG1002'}

    # Terpene synthase
    terpene_synthases = {'SMCOG1052'}

    # Cytochrome P450
    p450s = {'SMCOG1007', 'SMCOG1034'}

    # Hydrolase (a/b)
    hydrolases = {
        'SMCOG1036', 'SMCOG1066', 'SMCOG1115', 'SMCOG1179',
        'SMCOG1244'
    }

    # Reductase (FAD)
    reductases = {
        'SMCOG1001', 'SMCOG1161', 'SMCOG1083', 'SMCOG1196',
        'SMCOG1028', 'SMCOG1175', 'SMCOG1103', 'SMCOG1217',
        'SMCOG1039', 'SMCOG1176', 'SMCOG1159', 'SMCOG1219',
        'SMCOG1079', 'SMCOG1187', 'SMCOG1160'
    }

    # Dehydrogenase (SDR)
    dehydrogenases = {
        'SMCOG1001', 'SMCOG1117', 'SMCOG1040', 'SMCOG1222',
        'SMCOG1006', 'SMCOG1152', 'SMCOG1072', 'SMCOG1266',
        'SMCOG1017', 'SMCOG1181', 'SMCOG1095', 'SMCOG1293',
        'SMCOG1028', 'SMCOG1191', 'SMCOG1100'
    }

    # Dehalogenase (HAD)
    dehalogenases = {'SMCOG1119', 'SMCOG1262'}

    # Phosphatase
    phosphatases = {
        'SMCOG1012', 'SMCOG1145', 'SMCOG1101', 'SMCOG1249',
        'SMCOG1047', 'SMCOG1207', 'SMCOG1110', 'SMCOG1280',
        'SMCOG1060', 'SMCOG1219', 'SMCOG1123', 'SMCOG1283'
        'SMCOG1064', 'SMCOG1231'
    }

    # Transporter (MFS)
    transporters = {
        'SMCOG1000', 'SMCOG1044', 'SMCOG1085', 'SMCOG1096',
        'SMCOG1005', 'SMCOG1065', 'SMCOG1106', 'SMCOG1118',
        'SMCOG1020', 'SMCOG1067', 'SMCOG1131', 'SMCOG1166',
        'SMCOG1029', 'SMCOG1069', 'SMCOG1169', 'SMCOG1184',
        'SMCOG1035', 'SMCOG1074', 'SMCOG1202', 'SMCOG1205',
        'SMCOG1214', 'SMCOG1234', 'SMCOG1243', 'SMCOG1245',
        'SMCOG1252', 'SMCOG1254', 'SMCOG1288'
    }

    # Regulator
    regulators = {
        'SMCOG1008', 'SMCOG1136', 'SMCOG1071', 'SMCOG1224',
        'SMCOG1014', 'SMCOG1149', 'SMCOG1078', 'SMCOG1239',
        'SMCOG1015', 'SMCOG1167', 'SMCOG1112', 'SMCOG1260',
        'SMCOG1016', 'SMCOG1171', 'SMCOG1120', 'SMCOG1278',
        'SMCOG1031', 'SMCOG1195', 'SMCOG1125', 'SMCOG1284',
        'SMCOG1041', 'SMCOG1197', 'SMCOG1133', 'SMCOG1287',
        'SMCOG1057', 'SMCOG1201', 'SMCOG1135', 'SMCOG1058',
        'SMCOG1215'
    }

    # FAD-binding oxidases etc
    oxidases = {
        'SMCOG1138' 'SMCOG1055', 'SMCOG1138', 'SMCOG1217',
        'SMCOG1233', 'SMCOG1293'
    }

    # First, check PKS, NRPS and hybrids, since they require ALL domains
    if hybrid.issubset(domains):
        return 'PKS-NRPS'

    if pks.issubset(domains):
        return 'PKS'

    if nrps.issubset(domains):
        return 'NRPS'

    # Now, go through the rest just checking for any overlap and return
    # at the first match
    functions = [
        ('P450', p450s),
        ('Oxidase', oxidases),
        ('Regulator', regulators),
        ('Hydrolase', hydrolases),
        ('Reductase', reductases),
        ('Phosphatase', phosphatases),
        ('Transporter', transporters),
        ('Dehalogenase', dehalogenases),
        ('Dehydrogenase', dehydrogenases),
        ('Terpene synthase', terpene_synthases)
    ]

    for function_name, function_domains in functions:
        if not domains.isdisjoint(function_domains):
            return function_name

    # If nothing gets assigned, just say its hypothetical
    return 'Hypothetical'


def annotate_clusters(*clusters):
    """ Run hmmsearch on Gene objects in Clusters.

        1. Collect all Gene sequences in a NamedTemporaryFile
        2. Scan sequences for profiles in profiles.hmm
        3. Assign function to genes that meet thresholds

        Arguments:
            clusters (list): Cluster objects being scanned
    """
    annotations = defaultdict(dict)

    with NamedTemporaryFile() as fasta, \
            NamedTemporaryFile() as domtbl:

        # Build FASTA file, with headers corresponding to
        # indices of each Gene within its respective Cluster
        # for traceback purposes
        for i, cluster in enumerate(clusters):
            for j, gene in enumerate(cluster.genes):
                fasta.write(
                    f'>{i},{j}\n{gene.sequence}\n'.encode()
                )

        # Run hmmsearch on the FASTA tmpfile
        command = [
            HMMSEARCH_PATH,
            '--domtblout', domtbl.name,
            # '-E', '1E-16',  # antiSMASH setting
            HMMSEARCH_DATA,
            fasta.name
        ]

        subprocess.run(command, stdout=subprocess.DEVNULL)

        # Collect annotations for each gene
        for line in domtbl:
            line = line.decode()
            if line.startswith('#'):
                continue

            gene, *fields = line.split()
            # SMCOG1001:Description...
            smcog = fields[2].split(':')[0]

            domain = {
                smcog: (int(fields[18]),  # env from
                        int(fields[19]))  # env to
            }

            annotations[gene].update(domain)

    # Add functions to corresponding Gene objects
    for gene, domains in annotations.items():
        i, j = [int(x) for x in gene.split(',')]
        clusters[i].genes[j].function = assign_function(domains.keys())
        clusters[i].genes[j].smcog = domains
