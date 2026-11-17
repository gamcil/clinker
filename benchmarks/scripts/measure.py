#!/usr/bin/env python3
import argparse
import json
import subprocess
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--fixture", type=Path, required=True)
    parser.add_argument("--role", required=True)
    parser.add_argument("--clusters", type=int, required=True)
    parser.add_argument("--proteins", type=int, required=True)
    parser.add_argument("--mode", choices=("python", "rust", "rust-kmer"), required=True)
    parser.add_argument("--repeat", type=int, required=True)
    parser.add_argument("--threads", type=int, required=True)
    parser.add_argument("--identity", type=float, required=True)
    parser.add_argument("--kmer-size", type=int, required=True)
    parser.add_argument("--min-shared-kmers", type=int, required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    output_stem = Path(args.output)
    figure_path = Path(f"{output_stem}.figure.json")
    metrics_path = Path(f"{output_stem}.metrics.tmp")
    log_path = Path(f"{output_stem}.log")
    json_path = Path(f"{output_stem}.json")
    if args.mode == "python":
        command = ["/opt/python-clinker/bin/python", "/opt/benchmark/measure-python.py",
                   str(args.fixture), str(args.identity), str(args.threads),
                   str(figure_path), str(metrics_path)]
    else:
        command = ["/usr/local/bin/clinker-benchmark", str(args.fixture), args.mode,
                   str(args.identity), str(args.threads), str(args.kmer_size),
                   str(args.min_shared_kmers), str(figure_path), str(metrics_path)]
    completed = subprocess.run(command, capture_output=True)
    log_path.write_bytes(completed.stdout + completed.stderr)
    if completed.returncode:
        raise SystemExit(f"{args.mode} failed; see {log_path}")
    metrics = json.loads(metrics_path.read_text())
    metrics_path.unlink()

    record = {
        "dataset": args.fixture.name,
        "role": args.role,
        "clusters": args.clusters,
        "proteins": args.proteins,
        "mode": args.mode,
        "repeat": args.repeat,
        "identity": args.identity,
        "threads": args.threads,
        "kmer_size": args.kmer_size if args.mode == "rust-kmer" else None,
        "min_shared_kmers": args.min_shared_kmers if args.mode == "rust-kmer" else None,
        **metrics,
        "timing_boundary": "input-to-figure-json-v4",
        "alignments_performed": None,
    }
    json_path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
