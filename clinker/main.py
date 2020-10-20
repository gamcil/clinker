#!/usr/bin/env python3

"""
Command line interface

Cameron Gilchrist
"""

import argparse
import logging

from pathlib import Path

from clinker import align, plot
from clinker.classes import parse_files


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
LOG = logging.getLogger(__name__)


def clinker(
    files,
    identity=0.3,
    delimiter=None,
    decimals=2,
    output=None,
    force=False,
    hide_link_headers=False,
    hide_alignment_headers=False,
):
    """Entry point for running the script."""
    LOG.info("Starting clinker")

    # Check output file before doing anything else
    if output and Path(output).exists() and not force:
        LOG.error(f"File {output} already exists but --force not specified, exiting")
        return

    # Parse files, generate objects
    LOG.info("Parsing GenBank files")
    clusters = parse_files(files)

    # Align all clusters
    LOG.info("Aligning your clusters")
    globaligner = align.align_clusters(*clusters, cutoff=identity)

    # Generate results summary
    summary = globaligner.format(
        delimiter=delimiter,
        decimals=decimals,
        link_headers=not hide_link_headers,
        alignment_headers=not hide_alignment_headers,
    )
    if output:
        LOG.info(f"Writing alignments to {output}")
        with open(output, "w") as fp:
            fp.write(summary)
    else:
        print(summary)

    # Generate the SVG
    LOG.info("Plotting alignment...")
    plot.plot_clusters(globaligner)

    LOG.info("Done!")
    return globaligner


def get_parser():
    """Creates an ArgumentParser object."""
    parser = argparse.ArgumentParser(
        "clinker",
        description="clinker: Automatic creation of publication-ready"
        " gene cluster comparison figures.\n\n"
        "clinker generates gene cluster comparison figures from GenBank files."
        " It performs pairwise local or global alignments between every sequence"
        " in every unique pair of clusters and generates interactive, to-scale comparison figures"
        " using the clustermap.js library.",
        epilog="Example usage\n-------------\n"
        "Align clusters, plot results and print scores to screen:\n"
        "  $ clinker files/*.gbk\n\n"
        "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("files", help="Gene cluster GenBank files", nargs="+")

    alignment = parser.add_argument_group("Alignment options")
    alignment.add_argument(
        "-i",
        "--identity",
        help="Minimum alignment sequence identity",
        type=float,
        default=0.3
    )

    output = parser.add_argument_group("Output options")
    output.add_argument("-f", "--force", help="Overwrite previous output file", action="store_true")
    output.add_argument("-o", "--output", help="Save alignments to file")
    output.add_argument("-dl", "--delimiter", help="Character to delimit output by")
    output.add_argument("-dc", "--decimals", help="Number of decimal places in output", default=2)
    output.add_argument(
        "-hl",
        "--hide_link_headers",
        help="Hide alignment column headers",
        action="store_true",
    )
    output.add_argument(
        "-ha",
        "--hide_aln_headers",
        help="Hide alignment cluster name headers",
        action="store_true",
    )

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    clinker(
        args.files,
        identity=args.identity,
        delimiter=args.delimiter,
        decimals=args.decimals,
        output=args.output,
        force=args.force,
        hide_link_headers=args.hide_link_headers,
        hide_alignment_headers=args.hide_aln_headers,
    )


if __name__ == "__main__":
    main()
