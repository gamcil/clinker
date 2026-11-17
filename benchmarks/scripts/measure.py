#!/usr/bin/env python3
import argparse
import json
import subprocess
import time
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
    log_path = Path(f"{output_stem}.log")
    json_path = Path(f"{output_stem}.json")
    if args.mode == "python":
        command = ["/opt/python-clinker/bin/python", "/opt/benchmark/measure-python.py",
                   str(args.fixture), str(args.identity), str(args.threads)]
    else:
        inputs = sorted(args.fixture.glob("*.gbk"))
        if not inputs:
            raise SystemExit(f"No GenBank files found in {args.fixture}")
        command = ["/usr/local/bin/clinker-rust", *(str(path) for path in inputs),
                   "--plot-json", "--stats", "--identity", str(args.identity),
                   "--jobs", str(args.threads)]
        if args.mode == "rust-kmer":
            command += ["--prefilter", "--kmer-size", str(args.kmer_size),
                        "--min-shared-kmers", str(args.min_shared_kmers)]
    started = time.perf_counter()
    completed = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    elapsed = time.perf_counter() - started
    log_path.write_bytes(completed.stderr)
    if completed.returncode:
        raise SystemExit(f"{args.mode} failed; see {log_path}")
    figure_path.write_bytes(completed.stdout)
    figure = json.loads(completed.stdout)
    stats = next((json.loads(line.split("\t", 1)[1])
                  for line in completed.stderr.decode().splitlines()
                  if line.startswith("CLINKER_STATS\t")), None)
    metrics = {"elapsed_seconds": elapsed, "retained_links": len(figure["links"]),
               "groups": len(figure["groups"]), "json_bytes": len(completed.stdout),
               "alignments_performed": stats.get("alignments_performed") if stats else None}

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
        "timing_boundary": "cli-figure-json-v2",
    }
    json_path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
