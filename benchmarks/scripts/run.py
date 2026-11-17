#!/usr/bin/env python3
import argparse
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description='Run the Clinker benchmarks in this image.')
    parser.add_argument('--threads', type=int, default=4)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--fixture', action='append', help='Fixture ID; repeat to select multiple fixtures')
    parser.add_argument('--identity', type=float, default=0.3)
    parser.add_argument('--kmer-size', type=int, default=5)
    parser.add_argument('--min-shared-kmers', type=int, default=5)
    args = parser.parse_args()
    if not 1 <= args.threads <= 4 or args.repeats < 1:
        parser.error('threads must be 1–4 and repeats must be positive')
    workflow = Path('/opt/benchmark/workflow/Snakefile')
    fingerprint = Path('/opt/benchmark/build-fingerprint').read_text().strip()
    reference = os.environ.get('BENCHMARK_IMAGE', 'image digest not supplied')
    image = f'{reference}; build sha256:{fingerprint}'
    command = [
        'snakemake', '--snakefile', str(workflow), '--directory', '/work',
        '--cores', str(args.threads), '--resources', 'benchmark_slot=1',
        '--rerun-incomplete', '--config',
        f'threads={args.threads}', f'repeats={args.repeats}',
        f'fixtures={",".join(args.fixture) if args.fixture else "all"}',
        f'identity={args.identity}', f'kmer_size={args.kmer_size}',
        f'min_shared_kmers={args.min_shared_kmers}', f'image={image}',
    ]
    os.execvp(command[0], command)


if __name__ == '__main__':
    main()
