#!/usr/bin/env python3
import csv
import json
import statistics
import sys
from collections import defaultdict
from pathlib import Path


def median(records, field):
    values = [row[field] for row in records if row.get(field) is not None]
    return statistics.median(values) if values else ""


def main():
    *record_paths, output_path = map(Path, sys.argv[1:])
    records = [json.loads(path.read_text()) for path in record_paths]
    grouped = defaultdict(list)
    for row in records:
        grouped[(row["dataset"], row["mode"])].append(row)

    rows = []
    for dataset in sorted({row["dataset"] for row in records}):
        python_time = median(grouped[(dataset, "python")], "elapsed_seconds")
        rust_links = median(grouped[(dataset, "rust")], "retained_links")
        for mode in ("python", "rust", "rust-kmer"):
            values = grouped[(dataset, mode)]
            elapsed = median(values, "elapsed_seconds")
            links = median(values, "retained_links")
            rows.append({
                "dataset": dataset,
                "role": values[0]["role"],
                "clusters": values[0]["clusters"],
                "proteins": values[0]["proteins"],
                "mode": mode,
                "identity": values[0]["identity"],
                "threads": values[0]["threads"],
                "kmer_size": values[0]["kmer_size"],
                "min_shared_kmers": values[0]["min_shared_kmers"],
                "container": values[0]["container"],
                "python_clinker": values[0]["python_clinker"],
                "rust_clinker": values[0]["rust_clinker"],
                "python_version": values[0]["python_version"],
                "platform": values[0]["platform"],
                "architecture": values[0]["architecture"],
                "logical_cpus": values[0]["logical_cpus"],
                "median_elapsed_seconds": elapsed,
                "speedup_vs_python": python_time / elapsed,
                "median_retained_links": links,
                "links_vs_full_rust_percent": 100 * links / rust_links if rust_links else "",
                "median_potential_pairs": median(values, "potential_pairs"),
                "median_pairs_skipped_before_alignment": median(values, "pairs_skipped_before_alignment"),
                "median_alignments_performed": median(values, "alignments_performed"),
            })

    with output_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
