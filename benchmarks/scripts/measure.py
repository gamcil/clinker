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

    inputs = sorted(args.fixture.glob("*.gbk"))
    if not inputs:
        raise SystemExit(f"No GenBank files found in {args.fixture}")

    output_stem = Path(args.output)
    links_path = Path(f"{output_stem}.links")
    log_path = Path(f"{output_stem}.log")
    json_path = Path(f"{output_stem}.json")
    executable = "/opt/python-clinker/bin/clinker" if args.mode == "python" else "/usr/local/bin/clinker-rust"
    command = [
        executable,
        *(str(path) for path in inputs),
        "--identity", str(args.identity),
        "--jobs", str(args.threads),
        "--force",
    ]
    if args.mode == "python":
        command += ["--hide_link_headers", "--hide_aln_headers", "--delimiter", ","]
    else:
        command += ["--hide-link-headers", "--hide-alignment-headers"]
        if args.mode == "rust-kmer":
            command += [
                "--prefilter",
                "--kmer-size", str(args.kmer_size),
                "--min-shared-kmers", str(args.min_shared_kmers),
            ]
    command += ["--output", str(links_path)]

    started = time.perf_counter()
    completed = subprocess.run(command, text=True, capture_output=True)
    elapsed = time.perf_counter() - started
    log_path.write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise SystemExit(f"{args.mode} failed; see {log_path}")

    link_count = sum(bool(line.strip()) for line in links_path.read_text().splitlines())
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
        "elapsed_seconds": elapsed,
        "retained_links": link_count,
        "alignments_performed": None,
    }
    json_path.write_text(json.dumps(record, indent=2) + "\n")
    print(json.dumps(record), flush=True)


if __name__ == "__main__":
    main()
