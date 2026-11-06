#!/usr/bin/env bash
set -euo pipefail

script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
cd "$script_dir"

cores="${SNAKEMAKE_CORES:-${SLURM_CPUS_PER_TASK:-4}}"
exec snakemake \
  --use-conda \
  --conda-frontend conda \
  --cores "$cores" \
  --resources=benchmark_slot=1 \
  "$@"
