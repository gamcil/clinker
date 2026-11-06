#!/usr/bin/env python3
import argparse
import csv
import hashlib
import json
import os
import platform
import shutil
import statistics
import subprocess
import sys
import tarfile
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HERE = Path(__file__).resolve().parent
FIXTURES = {
    "bgcs_gcf_export_5": {"role": "small workload", "clusters": 15, "proteins": 221},
    "bgcs_gcf_export_6": {"role": "many clusters", "clusters": 150, "proteins": 1424},
    "bgcs_gcf_export_9": {"role": "many proteins", "clusters": 61, "proteins": 2110},
}
PYTHON_CLINKER_VERSION = "0.0.32"
IDENTITY = 0.3
THREADS = 4
REPEATS = 3
KMER_SIZE = 5
MIN_SHARED_KMERS = 8
MODES = ("python", "rust", "rust-kmer")


def run(command, **kwargs):
    return subprocess.run(command, check=True, text=True, **kwargs)


def prepare_tools():
    cache = HERE / ".cache"
    cache.mkdir(exist_ok=True)
    venv = cache / "python-venv"
    python_clinker = venv / "bin" / "clinker"
    if not python_clinker.exists():
        if venv.exists():
            shutil.rmtree(venv)
        python311 = shutil.which("python3.11")
        python = python311 or sys.executable
        run([python, "-m", "venv", str(venv)])
        run([
            str(venv / "bin" / "python"), "-m", "pip", "install",
            f"clinker=={PYTHON_CLINKER_VERSION}",
        ])

    cargo = shutil.which("cargo")
    if not cargo:
        raise SystemExit("Cargo is required to build the Rust CLI")
    run([
        cargo, "build", "--locked", "--release", "--manifest-path",
        str(ROOT / "Cargo.toml"), "-p", "clinker-cli",
    ])
    return python_clinker, ROOT / "target" / "release" / "clinker"


def prepare_fixture(name):
    archive = HERE / "fixtures" / f"{name}.tar.gz"
    expected = next(
        line.split()[0]
        for line in (HERE / "SHA256SUMS").read_text().splitlines()
        if line.endswith(f"fixtures/{name}.tar.gz")
    )
    digest = hashlib.sha256(archive.read_bytes()).hexdigest()
    if digest != expected:
        raise SystemExit(f"Fixture checksum failed: {archive}")

    fixture_dir = HERE / ".cache" / "data" / name
    files = sorted(fixture_dir.glob("*.gbk")) if fixture_dir.exists() else []
    if not files:
        data_dir = HERE / ".cache" / "data"
        data_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(archive, "r:gz") as fixture:
            fixture.extractall(data_dir)
        files = sorted(fixture_dir.glob("*.gbk"))
    if not files:
        raise SystemExit(f"No GenBank files found in fixture {name}")
    return files


def command_for(mode, executable, inputs, args, output):
    command = [str(executable), *(str(path) for path in inputs)]
    command += ["--identity", str(args.identity), "--jobs", str(args.threads), "--force"]
    if mode == "python":
        command += ["--hide_link_headers", "--hide_aln_headers"]
    else:
        command += ["--hide-link-headers", "--hide-alignment-headers", "--stats"]
        if mode == "rust-kmer":
            command += [
                "--prefilter", "--kmer-size", str(args.kmer_size),
                "--min-shared-kmers", str(args.min_shared_kmers),
            ]
    command += ["--output", str(output)]
    return command


def benchmark_one(mode, executable, inputs, args, dataset, repeat, results_dir):
    stem = f"{dataset}.{mode}.{repeat}"
    output = results_dir / "links" / f"{stem}.txt"
    log_path = results_dir / "logs" / f"{stem}.log"
    output.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    output.unlink(missing_ok=True)

    start = time.perf_counter()
    completed = subprocess.run(
        command_for(mode, executable, inputs, args, output),
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    elapsed = time.perf_counter() - start
    log_path.write_text(completed.stdout + completed.stderr)
    if completed.returncode:
        raise SystemExit(f"{mode} failed; see {log_path}")

    links = sum(bool(line.strip()) for line in output.read_text().splitlines())
    record = {
        "dataset": dataset,
        "role": FIXTURES[dataset]["role"],
        "clusters": FIXTURES[dataset]["clusters"],
        "proteins": FIXTURES[dataset]["proteins"],
        "mode": mode,
        "repeat": repeat,
        "elapsed_seconds": elapsed,
        "retained_links": links,
        "potential_pairs": None,
        "pairs_skipped_before_alignment": None,
        "alignments_performed": None,
    }
    for line in completed.stderr.splitlines():
        if line.startswith("CLINKER_STATS\t"):
            record.update(json.loads(line.split("\t", 1)[1]))
            break

    print(
        f"{dataset} {mode} run {repeat}: {elapsed:.3f}s, "
        f"{links} links" + (
            f", {record['alignments_performed']} alignments / "
            f"{record['potential_pairs']} potential pairs"
            if record["alignments_performed"] is not None else ""
        ),
        flush=True,
    )
    return record


def write_results(records, args, python_clinker, rust_clinker):
    results_dir = Path(args.output).resolve()
    results_dir.mkdir(parents=True, exist_ok=True)
    raw_fields = list(records[0])
    with (results_dir / "raw.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=raw_fields)
        writer.writeheader()
        writer.writerows(records)

    grouped = {}
    for record in records:
        grouped.setdefault((record["dataset"], record["mode"]), []).append(record)

    summary = []
    for dataset in args.datasets:
        python_time = statistics.median(
            row["elapsed_seconds"] for row in grouped[(dataset, "python")]
        )
        rust_links = statistics.median(
            row["retained_links"] for row in grouped[(dataset, "rust")]
        )
        for mode in MODES:
            rows = grouped[(dataset, mode)]
            elapsed = statistics.median(row["elapsed_seconds"] for row in rows)
            links = statistics.median(row["retained_links"] for row in rows)
            summary.append({
                "dataset": dataset,
                "role": FIXTURES[dataset]["role"],
                "clusters": FIXTURES[dataset]["clusters"],
                "proteins": FIXTURES[dataset]["proteins"],
                "mode": mode,
                "median_elapsed_seconds": elapsed,
                "median_retained_links": links,
                "speedup_vs_python": python_time / elapsed,
                "links_vs_full_rust_percent": 100 * links / rust_links,
                "median_potential_pairs": median(rows, "potential_pairs"),
                "median_pairs_skipped_before_alignment": median(
                    rows, "pairs_skipped_before_alignment"
                ),
                "median_alignments_performed": median(rows, "alignments_performed"),
            })

    with (results_dir / "summary.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(summary[0]))
        writer.writeheader()
        writer.writerows(summary)

    metadata = {
        "python_clinker": subprocess.check_output(
            [str(python_clinker), "--version"], text=True
        ).strip(),
        "rust_clinker": subprocess.check_output(
            [str(rust_clinker), "--version"], text=True
        ).strip(),
        "rustc": subprocess.check_output(["rustc", "--version"], text=True).strip(),
        "rust_source": subprocess.check_output(
            ["git", "-C", str(ROOT), "describe", "--always", "--dirty"], text=True
        ).strip(),
        "identity_cutoff": args.identity,
        "threads": args.threads,
        "repeats": args.repeats,
        "kmer_size": args.kmer_size,
        "min_shared_kmers": args.min_shared_kmers,
        "platform": platform.platform(),
        "processor": platform.processor(),
        "logical_cpus": os.cpu_count(),
        "datasets": {dataset: FIXTURES[dataset] for dataset in args.datasets},
    }
    (results_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def median(rows, key):
    values = [row[key] for row in rows if row[key] is not None]
    return statistics.median(values) if values else None


def parse_args():
    parser = argparse.ArgumentParser(description="Run Clinker runtime benchmarks")
    parser.add_argument("--dataset", choices=("all", *FIXTURES), default="all")
    parser.add_argument("--threads", type=int, default=THREADS)
    parser.add_argument("--repeats", type=int, default=REPEATS)
    parser.add_argument("--identity", type=float, default=IDENTITY)
    parser.add_argument("--kmer-size", type=int, default=KMER_SIZE)
    parser.add_argument("--min-shared-kmers", type=int, default=MIN_SHARED_KMERS)
    parser.add_argument("--output", default=str(HERE / "results"))
    args = parser.parse_args()
    args.datasets = list(FIXTURES) if args.dataset == "all" else [args.dataset]
    return args


def main():
    args = parse_args()
    python_clinker, rust_clinker = prepare_tools()
    records = []
    for dataset in args.datasets:
        inputs = prepare_fixture(dataset)
        for mode in MODES:
            executable = python_clinker if mode == "python" else rust_clinker
            for repeat in range(1, args.repeats + 1):
                records.append(
                    benchmark_one(mode, executable, inputs, args, dataset, repeat, Path(args.output))
                )
    write_results(records, args, python_clinker, rust_clinker)
    print(f"Results written to {Path(args.output).resolve()}")


if __name__ == "__main__":
    main()
