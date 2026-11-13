#!/usr/bin/env python3
import json
import logging
import sys
import time
from pathlib import Path

from clinker.align import align_clusters
from clinker.classes import parse_files

def main():
    logging.getLogger("clinker").setLevel(logging.ERROR)
    fixture, identity, threads, figure_path, metrics_path = sys.argv[1:]
    inputs = sorted(str(path) for path in Path(fixture).glob('*.gbk'))
    if not inputs:
        raise SystemExit(f'No GenBank files found in {fixture}')

    started = time.perf_counter()
    clusters = parse_files(inputs)
    analysis = align_clusters(*clusters, cutoff=float(identity), jobs=int(threads))
    analysis.build_gene_groups()
    figure = analysis.to_data()
    data = json.dumps(figure, separators=(',', ':')).encode('utf-8')
    elapsed = time.perf_counter() - started

    Path(figure_path).write_bytes(data)
    Path(metrics_path).write_text(json.dumps({
        'elapsed_seconds': elapsed, 'retained_links': len(figure['links']),
        'groups': len(figure['groups']), 'json_bytes': len(data),
    }))


if __name__ == "__main__":
    main()
