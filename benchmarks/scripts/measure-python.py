#!/usr/bin/env python3
import json
import logging
import sys
from pathlib import Path

from clinker.align import align_clusters
from clinker.classes import parse_files

def main():
    logging.getLogger("clinker").setLevel(logging.ERROR)
    fixture, identity, threads, output = sys.argv[1:]
    inputs = sorted(str(path) for path in Path(fixture).glob('*.gbk'))
    if not inputs:
        raise SystemExit(f'No GenBank files found in {fixture}')

    clusters = parse_files(inputs)
    analysis = align_clusters(*clusters, cutoff=float(identity), jobs=int(threads))
    analysis.build_gene_groups()
    if analysis.alignments:
        print(analysis.format())
    figure = analysis.to_data()
    Path(output).write_text(json.dumps(figure, indent=2))


if __name__ == "__main__":
    main()
