#!/usr/bin/env python3

"""
Command line interface

Cameron Gilchrist
"""

import argparse
import logging

from pathlib import Path

from clinker import align
from clinker.plot import plot_clusters
from clinker.classes import find_files, parse_files


logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s] %(levelname)s - %(message)s",
    datefmt="%H:%M:%S"
)
LOG = logging.getLogger(__name__)


def clinker(
    files,
    session=None,
    identity=0.3,
    delimiter=None,
    decimals=2,
    plot=None,
    output=None,
    force=False,
    no_align=False,
    hide_link_headers=False,
    hide_alignment_headers=False,
    use_file_order=False,
    json_indent=None,
    jobs=None,
):
    """Entry point for running the script."""
    LOG.info("Starting clinker")

    load_session = session and Path(session).exists()

    if load_session:
        LOG.info("Loading session from: %s", session)
        with open(session) as fp:
            globaligner = align.Globaligner.from_json(fp)
        if files:
            paths = find_files(files)
            if not paths:
                LOG.error("No files found")
                raise SystemExit
            LOG.info("Parsing GenBank files: %s", paths)
            clusters = parse_files(paths)

            LOG.info("Adding clusters to loaded session and aligning")
            globaligner.add_clusters(*clusters)
            globaligner.align_stored_clusters(cutoff=identity, jobs=jobs)
            load_session = False
    else:
        # Parse files, generate objects
        paths = find_files(files)
        if not paths:
            LOG.error("No files found")
            raise SystemExit
        LOG.info("Parsing GenBank files: %s", paths)
        clusters = parse_files(paths)

        # Align all clusters
        if no_align:
            globaligner = align.Globaligner()
            globaligner.add_clusters(*clusters)
        elif len(clusters) == 1:
            globaligner = align.align_clusters(clusters[0], jobs=1)
        else:
            LOG.info("Starting cluster alignments")
            globaligner = align.align_clusters(*clusters, cutoff=identity, jobs=jobs)

    if globaligner.alignments:
        LOG.info("Generating results summary...")
        summary = globaligner.format(
            delimiter=delimiter,
            decimals=decimals,
            link_headers=not hide_link_headers,
            alignment_headers=not hide_alignment_headers,
        )
        if output:
            if (output and Path(output).exists() and not force):
                print(summary)
                LOG.warn("File %s already exists but --force was not specified", output)
            else:
                LOG.info("Writing alignments to: %s", output)
                with open(output, "w") as fp:
                    fp.write(summary)
        else:
            print(summary)
    else:
        LOG.info("No alignments were generated")

    if session and not load_session:
        LOG.info("Saving session to: %s", session)
        with open(session, "w") as fp:
            globaligner.to_json(fp, indent=json_indent)

    # Generate the SVG
    if plot:
        LOG.info("Building clustermap.js visualisation")
        if isinstance(plot, str):
            LOG.info("Writing to: %s", plot)
        plot_clusters(
            globaligner,
            output=None if plot is True else plot,
            use_file_order=use_file_order,
        )

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
        "Only save gene-gene links when identity is over 50%:\n"
        "  $ clinker files/*.gbk -i 0.5\n\n"
        "Save an alignment session for later:\n"
        "  $ clinker files/*.gbk -s session.json\n\n"
        "Save alignments to file, in comma-delimited format, with 4 decimal places:\n"
        "  $ clinker files/*.gbk -o alignments.csv -dl \",\" -dc 4\n\n"
        "Generate visualisation:\n"
        "  $ clinker files/*.gbk -p\n\n"
        "Save visualisation as a static HTML document:\n"
        "  $ clinker files/*.gbk -p plot.html\n\n"
        "Cameron Gilchrist, 2020",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("files", help="Gene cluster GenBank files", nargs="+")

    alignment = parser.add_argument_group("Alignment options")
    alignment.add_argument(
        "-na",
        "--no_align",
        help="Do not align clusters",
        action="store_true",
    )
    alignment.add_argument(
        "-i",
        "--identity",
        help="Minimum alignment sequence identity [default: 0.3]",
        type=float,
        default=0.3
    )
    alignment.add_argument(
        "-j",
        "--jobs",
        help="Number of alignments to run in parallel (0 to use the number of CPUs) [default: 0]",
        type=int,
        default=0,
    )

    output = parser.add_argument_group("Output options")
    output.add_argument("-s", "--session", help="Path to clinker session")
    output.add_argument("-ji", "--json_indent", type=int, help="Number of spaces to indent JSON [default: none]")
    output.add_argument("-f", "--force", help="Overwrite previous output file", action="store_true")
    output.add_argument("-o", "--output", help="Save alignments to file")
    output.add_argument(
        "-p",
        "--plot",
        nargs="?",
        const=True,
        default=False,
        help="Plot cluster alignments using clustermap.js. If a path is given,"
        " clinker will generate a portable HTML file at that path. Otherwise,"
        " the plot will be served dynamically using Python's HTTP server."
    )
    output.add_argument("-dl", "--delimiter", help="Character to delimit output by [default: human readable]")
    output.add_argument("-dc", "--decimals", help="Number of decimal places in output [default: 2]", default=2)
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

    viz = parser.add_argument_group("Visualisation options")
    viz.add_argument(
        "-ufo",
        "--use_file_order",
        action="store_true",
        help="Display clusters in order of input files"
    )

    return parser


def main():
    parser = get_parser()
    args = parser.parse_args()
    clinker(
        args.files,
        session=args.session,
        json_indent=args.json_indent,
        identity=args.identity,
        delimiter=args.delimiter,
        decimals=args.decimals,
        plot=args.plot,
        output=args.output,
        force=args.force,
        no_align=args.no_align,
        hide_link_headers=args.hide_link_headers,
        hide_alignment_headers=args.hide_aln_headers,
        use_file_order=args.use_file_order,
        jobs=args.jobs if args.jobs > 0 else None
    )


if __name__ == "__main__":
    main()
