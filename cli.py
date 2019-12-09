#!/usr/bin/env python3

"""
Command line interface

Cameron Gilchrist
"""

import argparse

from Bio import SeqIO

import globaligner
import visualise

from classes import Cluster


def run(args):
    """ Entry point for running the script.
    """
    # Parse files, generate objects
    clusters = []
    print("Parsing GenBank files")
    for file in args.files:
        records = [record for record in SeqIO.parse(file, "genbank")]
        cluster = Cluster.from_seqrecords(*records, mode=args.mode)
        clusters.append(cluster)

    # Scan genes for smCOGs if protein mode specified
    if args.mode == "protein":
        print("Scanning for biosynthetic function")
        import annotate

        # annotate.find_biosynthetic_functions(clusters)

    # Align all clusters
    print("Aligning your clusters")
    aligner = globaligner.align_clusters(*clusters)

    # Generate the SVG
    print("Generating the figure")
    figure = visualise.Figure(aligner=aligner, width=args.width).render("maximize")

    # Write output files
    with open(args.output, "w") as svg:

        # SVG is just one big string
        svg.write(figure)

    print(aligner)
    print("Done!")


def main():
    """ Collect args, then launch script with run()
    """
    parser = argparse.ArgumentParser(
        prog="crosslinker",
        description="Automatic creation of high quality comparison"
        " figures for biosynthetic gene clusters (BGCs)."
        " Given GenBank files containing BGCs, this"
        " package will parse them for genes/proteins,"
        " scan for biosynthetic functions (protein only),"
        " perform pairwise alignments of all sequences,"
        " and generate a to-scale comparison figure.",
        epilog="Cameron Gilchrist 2019",
    )

    # Required args
    required = parser.add_argument_group("Required arguments")
    required.add_argument(
        "mode", choices=["protein", "gene"], help="Running mode of the script"
    )
    required.add_argument(
        "output",
        help="Base name for output files (e.g. foo -> foo.svg, foo.homologies)",
    )
    required.add_argument(
        "--files", help="Gene cluster GenBank files", nargs="+", required=True
    )

    # Optional args
    optional = parser.add_argument_group("Optional arguments")
    optional.add_argument(
        "--colours", help="Custom colour scheme for biosynthetic functions"
    )
    optional.add_argument(
        "--width", help="Width (pixels) of the SVG figure", default=1000
    )
    optional.add_argument(
        "--spacing",
        help="Number of vertical space (pixels) between each cluster",
        default=40,
    )
    optional.add_argument(
        "--ordered",
        help="Place clusters in the figure in the same order as your"
        " input files. By default, the first file will always"
        " appear first, but subsequent clusters will be shown"
        " in homology score order (detailed in README)",
        action="store_true",
    )

    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
