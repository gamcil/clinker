#!/usr/bin/env python3
import csv
import json
import os
import platform
import statistics
import subprocess
import sys
from collections import defaultdict
from pathlib import Path


def write_csv(path, rows):
    with Path(path).open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    image, output_dir, *paths = sys.argv[1:]
    output_dir = Path(output_dir)
    records = [json.loads(Path(path).read_text()) for path in paths]
    browser_versions = sorted({row.pop('browser_version') for row in records if row.get('browser_version')})
    metadata = {
        'container': image,
        'python_clinker': subprocess.check_output(['/opt/python-clinker/bin/clinker', '--version'], text=True).strip(),
        'rust_clinker': subprocess.check_output(['/usr/local/bin/clinker-rust', '--version'], text=True).strip(),
        'web_build': json.loads(Path('/opt/clinker-web/build.json').read_text())['id'],
        'browser_versions': browser_versions,
        'python_version': platform.python_version(),
        'platform': platform.platform(),
        'architecture': platform.machine(),
        'logical_cpus': os.cpu_count(),
    }
    (output_dir / 'run-metadata.json').write_text(json.dumps(metadata, indent=2) + '\n')

    fields = ('dataset', 'role', 'clusters', 'proteins', 'mode', 'repeat', 'identity',
              'threads', 'kmer_size', 'min_shared_kmers', 'elapsed_seconds',
              'retained_links', 'alignments_performed', 'timing_boundary')
    records = [{field: row.get(field) for field in fields} for row in records]
    records.sort(key=lambda row: (row['dataset'], row['mode'], row['repeat']))
    write_csv(output_dir / 'runs.csv', records)

    grouped = defaultdict(list)
    for row in records:
        grouped[row['dataset'], row['mode']].append(row)

    summaries = []
    for dataset in sorted({row['dataset'] for row in records}):
        python_time = statistics.median(row['elapsed_seconds'] for row in grouped[dataset, 'python'])
        rust_links = statistics.median(row['retained_links'] for row in grouped[dataset, 'rust'])
        for mode in ('python', 'rust', 'rust-kmer', 'wasm', 'wasm-kmer'):
            group = grouped[dataset, mode]
            first = group[0]
            elapsed = statistics.median(row['elapsed_seconds'] for row in group)
            links = statistics.median(row['retained_links'] for row in group)
            summaries.append({
                **{key: first[key] for key in ('dataset', 'role', 'clusters', 'proteins', 'mode',
                                                  'identity', 'threads', 'kmer_size', 'min_shared_kmers')},
                'repeats': len(group),
                'median_elapsed_seconds': elapsed,
                'speedup_vs_python': python_time / elapsed,
                'median_retained_links': links,
                'links_vs_full_rust_percent': 100 * links / rust_links if rust_links else None,
            })
    write_csv(output_dir / 'summary.csv', summaries)


if __name__ == '__main__':
    main()
