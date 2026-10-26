#!/usr/bin/env python3
"""Compare the legacy Python, native Rust and browser WASM analysis paths."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import select
import signal
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
NATIVE = ROOT / "target/release/examples/benchmark"
FIELDS = ["dataset", "configuration", "iteration", "workers", "identity", "status",
          "kmer_size", "min_shared_kmers",
          "parse_ms", "prepare_ms", "prefilter_ms", "alignment_ms", "total_ms",
          "clusters", "proteins", "amino_acids", "possible_pairs", "eligible_pairs",
          "candidate_pairs", "aligned_pairs", "retained_links", "error"]
CONFIGS = ["python", "native", "wasm", "native-prefilter", "wasm-prefilter"]


def sha(data):
    return hashlib.sha256(data).hexdigest()


def stop(process):
    if process.poll() is None:
        os.killpg(process.pid, signal.SIGTERM)
        try:
            process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGKILL)
    process.wait()


def execute(command, request, timeout, log):
    with log.open("w") as stderr:
        process = subprocess.Popen(command, cwd=ROOT, stdin=subprocess.PIPE,
                                   stdout=subprocess.PIPE, stderr=stderr, text=True,
                                   start_new_session=True)
        try:
            output, _ = process.communicate(json.dumps(request), timeout=timeout)
        except subprocess.TimeoutExpired:
            stop(process)
            raise RuntimeError("timeout after {} seconds; see {}".format(timeout, log))
        except BaseException:
            stop(process)
            raise
    if process.returncode:
        raise RuntimeError("process exited {}; see {}".format(process.returncode, log))
    return json.loads(output)


class Browser:
    def __init__(self, output, timeout):
        self.timeout = timeout
        self.log = (output / "browser.stderr.log").open("w")
        self.process = subprocess.Popen(["node", str(ROOT / "bench/browser.mjs")], cwd=ROOT,
                                        stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                        stderr=self.log, text=True, bufsize=1,
                                        start_new_session=True)
        try:
            self.metadata = self.read()
            if not self.metadata.get("ready"):
                raise RuntimeError("browser did not start")
        except Exception:
            self.close()
            raise

    def read(self):
        if not select.select([self.process.stdout], [], [], self.timeout)[0]:
            raise RuntimeError("browser timeout; see browser.stderr.log")
        line = self.process.stdout.readline()
        if not line:
            raise RuntimeError("browser exited; see browser.stderr.log")
        result = json.loads(line)
        if result.get("error"):
            raise RuntimeError(result["error"])
        return result

    def run(self, request):
        self.process.stdin.write(json.dumps(request) + "\n")
        self.process.stdin.flush()
        return self.read()

    def close(self):
        stop(self.process)
        self.log.close()


def canonical(result):
    # Ignore within-file feature order (e.g. circular-origin rotation), but
    # preserve cluster/file membership and duplicated sequences.
    proteins = [sorted(cluster) for cluster in result["proteins"]]
    if any(not protein for cluster in proteins for protein in cluster):
        raise RuntimeError("empty protein translation; resolve parsing before benchmarking")
    return proteins


def version(command):
    return subprocess.check_output(command, cwd=ROOT, text=True).strip()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--datasets", nargs="+", default=["bgcs", "bgcs_gcf_export_1"])
    parser.add_argument("--repeats", type=int, default=5)
    parser.add_argument("--workers", type=int, default=4)
    parser.add_argument("--identity", type=float, default=0.3)
    parser.add_argument("--configurations", nargs="+", choices=CONFIGS, default=CONFIGS)
    parser.add_argument("--kmer-size", type=int, default=3)
    parser.add_argument("--min-shared-kmers", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--timeout", type=float, default=3600)
    parser.add_argument("--python", default=sys.executable, help="Interpreter with v1 dependencies")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--legacy-source", type=Path, required=True, help="Python maintenance checkout or extracted release")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--prepare-inputs", action="store_true", help="Create separate copies excluding documented clipped boundary CDS fragments")
    parser.add_argument("--smoke", action="store_true", help="Two smallest files per dataset; one repetition")
    args = parser.parse_args()
    args.legacy_source = args.legacy_source.resolve()
    if not (args.legacy_source / "clinker/align.py").is_file():
        parser.error("--legacy-source must point to a Python maintenance checkout or extracted release")
    if min(args.repeats, args.workers, args.timeout) <= 0 or not 0 < args.identity <= 1:
        parser.error("positive repeats/workers/timeout and identity in (0,1] required")
    if not 1 <= args.kmer_size <= 8 or args.min_shared_kmers < 1:
        parser.error("kmer-size must be in [1,8] and min-shared-kmers must be positive")
    if len(set(args.configurations)) != len(args.configurations):
        parser.error("configurations must be unique")
    if args.smoke:
        args.repeats = 1
    if not NATIVE.is_file() or not (ROOT / "web/pkg/clinker_wasm_bg.wasm").is_file():
        parser.error("build the native example and web WASM first; see bench/README.md")
    output = (args.output or ROOT / "benchmark-results" / time.strftime("%Y%m%d-%H%M%S")).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        parser.error("output directory must be empty")
    if args.prepare_inputs:
        sources = [(ROOT / directory).resolve() for directory in args.datasets]
        source_manifest = {str(source): [{"name": p.name, "sha256": sha(p.read_bytes())}
                                         for p in sorted(source.iterdir()) if p.suffix.lower() in {".gbk", ".gb", ".gbff", ".genbank"}]
                           for source in sources}
        (output / "source-datasets.json").write_text(json.dumps(source_manifest, indent=2) + "\n")
        prepared = execute([args.python, str(ROOT / "bench/prepare_inputs.py")],
                           {"datasets": list(map(str, sources)), "output": str(output / "inputs")},
                           args.timeout, output / "prepare.stderr.log")
        (output / "preparation.json").write_text(json.dumps(prepared, indent=2) + "\n")
        args.datasets = prepared["datasets"]
        print("Prepared copies; removed {} boundary CDS fragments. Originals unchanged.".format(len(prepared["removed"])), flush=True)
    datasets = {}
    manifest = {}
    for directory in args.datasets:
        directory = (ROOT / directory).resolve()
        files = sorted(p for p in directory.iterdir() if p.suffix.lower() in {".gbk", ".gb", ".gbff", ".genbank"})
        if args.smoke:
            files = sorted(sorted(files, key=lambda p: (p.stat().st_size, p.name))[:2])
        if len(files) < 2 or directory.name in datasets:
            parser.error("each dataset needs a unique directory name and at least two GenBank files")
        datasets[directory.name] = [str(p) for p in files]
        manifest[directory.name] = [{"path": str(p), "bytes": p.stat().st_size,
                                     "sha256": sha(p.read_bytes())} for p in files]
    (output / "datasets.json").write_text(json.dumps(manifest, indent=2) + "\n")
    env = os.environ.copy()
    for name in ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"]:
        os.environ[name] = "1"
    browser = None
    try:
        browser = Browser(output, args.timeout)
        metadata = {"arguments": vars(args), "platform": platform.platform(), "cpu": platform.processor(),
                    "logical_cpus": os.cpu_count(), "browser": browser.metadata,
                    "git_commit": version(["git", "rev-parse", "HEAD"]),
                    "git_status": version(["git", "status", "--short"]),
                    "python": version([args.python, "--version"]),
                    "python_packages": version([args.python, "-m", "pip", "freeze"]),
                    "rustc": version(["rustc", "-Vv"]),
                    "native_sha256": sha(NATIVE.read_bytes()),
                    "wasm_sha256": sha((ROOT / "web/pkg/clinker_wasm_bg.wasm").read_bytes()),
                    "adapter_sha256": {p: sha((ROOT / p).read_bytes()) for p in
                        ["bench/run.py", "bench/python_v1.py", "bench/browser.mjs", "bench/wasm.mjs",
                         "bench/prepare_inputs.py", "crates/clinker-core/examples/benchmark.rs",
                         "web/alignment.mjs", "web/worker.js", "Cargo.lock"]},
                    "legacy_source_sha256": {p: sha((args.legacy_source / p).read_bytes()) for p in
                        ["clinker/align.py", "clinker/classes.py", "clinker/formatters.py"]},
                    "cache_policy": "warm filesystem; fresh analysis and alignment pool each run",
                    "prefilter": {"k": args.kmer_size, "count": args.min_shared_kmers}}
        (output / "metadata.json").write_text(json.dumps(metadata, indent=2, default=str) + "\n")

        def run(config, files, inspect=False, label="run"):
            request = {"legacy_source": str(args.legacy_source), "files": files, "workers": args.workers, "identity": args.identity,
                       "prefilter": config.endswith("-prefilter"), "inspect": inspect,
                       "kmer_size": args.kmer_size, "min_shared_kmers": args.min_shared_kmers}
            if config.startswith("wasm"):
                return browser.run(request)
            command = [args.python, str(ROOT / "bench/python_v1.py")] if config == "python" else [str(NATIVE)]
            return execute(command, request, args.timeout, output / (label + ".stderr.log"))

        validation = {}
        for dataset, files in datasets.items():
            print("Validating {} ({} files)...".format(dataset, len(files)), flush=True)
            outputs = {config: canonical(run(config, files, True, "validate-" + dataset + "-" + config))
                       for config in ["python", "native", "wasm"]}
            (output / ("proteins-" + dataset + ".json")).write_text(json.dumps(outputs) + "\n")
            if outputs["python"] != outputs["native"] or outputs["native"] != outputs["wasm"]:
                raise RuntimeError("{}: parsed proteins differ; see proteins-{}.json. No timed runs performed.".format(dataset, dataset))
            counts = list(map(len, outputs["native"]))
            validation[dataset] = {"clusters": len(counts), "proteins": sum(counts),
                                   "amino_acids": sum(len(p) for cluster in outputs["native"] for p in cluster),
                                   "possible_pairs": sum(a*b for i, a in enumerate(counts) for b in counts[i+1:]),
                                   "sha256": sha(json.dumps(outputs["native"]).encode())}
        (output / "validation.json").write_text(json.dumps(validation, indent=2) + "\n")
        if args.validate_only:
            print("All protein inputs match. Results: {}".format(output))
            return
        rng = random.Random(args.seed)
        schedule = []
        for iteration in range(1, args.repeats + 1):
            block = [(dataset, config, iteration) for dataset in datasets for config in args.configurations]
            rng.shuffle(block)
            schedule.extend(block)
        (output / "schedule.json").write_text(json.dumps(schedule, indent=2) + "\n")
        rows = []
        with (output / "runs.csv").open("w", newline="") as handle:
            writer = csv.DictWriter(handle, FIELDS)
            writer.writeheader()
            handle.flush()
            for dataset, config, iteration in schedule:
                label = "{}-{}-{}".format(dataset, config, iteration)
                print("Running " + label, flush=True)
                row = {"dataset": dataset, "configuration": config, "iteration": iteration,
                       "workers": args.workers, "identity": args.identity, "status": "ok",
                       "kmer_size": args.kmer_size if config.endswith("-prefilter") else None,
                       "min_shared_kmers": args.min_shared_kmers if config.endswith("-prefilter") else None}
                try:
                    row.update(run(config, datasets[dataset], label=label))
                    if row["proteins"] != validation[dataset]["proteins"]:
                        raise RuntimeError("protein count changed after validation")
                except Exception as error:
                    row.update(status="error", error=str(error))
                writer.writerow(row)
                handle.flush()
                rows.append(row)
                if row["status"] != "ok":
                    raise RuntimeError("{} failed: {} (partial CSV saved)".format(label, row["error"]))
                print("  {:.3f} s; {} alignments; {} retained links".format(
                    row["total_ms"] / 1000, row["aligned_pairs"], row["retained_links"]), flush=True)
        for dataset in datasets:
            for iteration in range(1, args.repeats + 1):
                block = {r["configuration"]: r for r in rows if r["dataset"] == dataset and r["iteration"] == iteration}
                if {"native", "wasm"} <= block.keys() and block["native"]["aligned_pairs"] != block["wasm"]["aligned_pairs"]:
                    raise RuntimeError(dataset + ": native/WASM unfiltered pair counts differ")
                if {"native-prefilter", "wasm-prefilter"} <= block.keys() and block["native-prefilter"]["candidate_pairs"] != block["wasm-prefilter"]["candidate_pairs"]:
                    raise RuntimeError(dataset + ": native/WASM prefilter candidate counts differ")
        # Catch file edits during a long benchmark rather than accepting mixed inputs.
        for files in manifest.values():
            for item in files:
                if sha(Path(item["path"]).read_bytes()) != item["sha256"]:
                    raise RuntimeError("input changed during benchmark: " + item["path"])
        with (output / "summary.csv").open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["dataset", "configuration", "kmer_size", "min_shared_kmers", "runs", "median_total_ms", "min_total_ms", "max_total_ms"])
            for dataset in datasets:
                for config in args.configurations:
                    times = [r["total_ms"] for r in rows if r["dataset"] == dataset and r["configuration"] == config]
                    writer.writerow([dataset, config,
                                     args.kmer_size if config.endswith("-prefilter") else None,
                                     args.min_shared_kmers if config.endswith("-prefilter") else None,
                                     len(times), statistics.median(times), min(times), max(times)])
        print("Results: {}".format(output))
    finally:
        if browser:
            browser.close()
        os.environ.clear()
        os.environ.update(env)


if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print("Benchmark cancelled; partial results retained.", file=sys.stderr)
        sys.exit(130)
    except (RuntimeError, OSError, ValueError) as error:
        print("Benchmark failed: " + str(error), file=sys.stderr)
        sys.exit(1)
