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


def figure_members(path):
    figure = json.loads(Path(path).read_text())
    genes = {}
    keys = set()
    for cluster in figure['clusters']:
        for locus in cluster['loci']:
            for gene in locus['genes']:
                key = (cluster['name'], locus['name'], gene['label'])
                if key in keys:
                    raise ValueError('gene labels are not unique within a locus')
                keys.add(key)
                genes[gene['uid']] = key
    links = {tuple(sorted((genes[link['query']['uid']], genes[link['target']['uid']])))
             for link in figure['links']}
    groups = {tuple(sorted(genes[uid] for uid in group['genes'])) for group in figure['groups']}
    return links, groups, keys


def main():
    image, output_dir, *paths = sys.argv[1:]
    output_dir = Path(output_dir)
    records = [json.loads(Path(path).read_text()) for path in paths]
    by_run = defaultdict(dict)
    for row, path in zip(records, paths):
        figure_path = str(Path(path).with_suffix('')) + '.figure.json'
        by_run[row['dataset'], row['repeat']][row['mode']] = (row, figure_path)
    for modes in by_run.values():
        try:
            reference_links, reference_groups, reference_genes = figure_members(modes['rust'][1])
        except (KeyError, ValueError):
            for row, _ in modes.values():
                row['comparison_note'] = 'Rust reference could not be compared'
            continue
        for row, path in modes.values():
            try:
                links, groups, genes = figure_members(path)
                row['genes_only_in_run'] = len(genes - reference_genes)
                row['genes_only_in_rust'] = len(reference_genes - genes)
                if genes != reference_genes:
                    row['comparison_note'] = 'Gene identities differ; shared counts may be incomplete'
                row['links_shared_with_rust'] = len(links & reference_links)
                row['groups_shared_with_rust'] = len(groups & reference_groups)
            except (KeyError, ValueError) as error:
                row['comparison_note'] = str(error)
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
              'retained_links', 'groups', 'json_bytes', 'links_shared_with_rust',
              'groups_shared_with_rust', 'genes_only_in_run', 'genes_only_in_rust', 'comparison_note',
              'alignments_performed', 'timing_boundary')
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
                'median_groups': statistics.median(row['groups'] for row in group),
                'median_links_shared_with_rust': (statistics.median(row['links_shared_with_rust'] for row in group)
                                                  if all(row['links_shared_with_rust'] is not None for row in group) else None),
                'median_groups_shared_with_rust': (statistics.median(row['groups_shared_with_rust'] for row in group)
                                                   if all(row['groups_shared_with_rust'] is not None for row in group) else None),
                'links_vs_full_rust_percent': 100 * links / rust_links if rust_links else None,
            })
    write_csv(output_dir / 'summary.csv', summaries)


if __name__ == '__main__':
    main()
