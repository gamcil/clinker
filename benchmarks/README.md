# Clinker runtime benchmarks

This standalone Nextflow workflow compares pinned published releases of Python
Clinker and the Rust Clinker CLI. It runs every fixture with the Python CLI,
the Rust CLI, and the Rust CLI with the k-mer prefilter enabled. Each run is an
independent cached task; the timed interval covers the CLI invocation only.

## Run

Requirements: Nextflow 25.10.2 or newer, Java, and Docker. For an HPC system,
use Apptainer instead of Docker.

```sh
docker build -t ghcr.io/gamcil/clinker-benchmarks:0.1.0 .
nextflow run main.nf -profile docker -resume
```

To use Apptainer on a cluster, publish the built image to GHCR first:

```sh
docker login ghcr.io
docker push ghcr.io/gamcil/clinker-benchmarks:0.1.0
nextflow run main.nf -profile apptainer -resume
```

For a release benchmark, inspect the pushed image and pass its immutable digest
with `--container ghcr.io/gamcil/clinker-benchmarks@sha256:<digest>`.

Keep Nextflow's `work/` and `.nextflow/` directories between attempts so
`-resume` can reuse completed fixture extractions and benchmark runs. Results,
logs, link files, a trace, execution report, and timeline are written under
`results/`.

To select one fixture or alter repetitions and threads:

```sh
nextflow run main.nf -profile docker -resume \
  --fixtures 'fixtures/bgcs_gcf_export_5.tar.gz' \
  --repeats 1 --threads 2
```

## Fixtures

Five fixed GenBank fixture archives are in `fixtures/`. Their checksums are in
`SHA256SUMS`. Two are fixed samples from the larger export fixtures: export 10
contains 50 clusters sampled from export 6; export 11 combines 50 clusters
from export 6 with 25 from export 9.

| Fixture | Workload | Clusters | Proteins |
| --- | --- | ---: | ---: |
| `bgcs_gcf_export_5` | Small | 15 | 221 |
| `bgcs_gcf_export_10` | Medium, balanced | 50 | 490 |
| `bgcs_gcf_export_9` | Protein-heavy | 61 | 2,110 |
| `bgcs_gcf_export_11` | Mixed | 75 | 1,359 |
| `bgcs_gcf_export_6` | Cluster-heavy | 150 | 1,424 |

`fixtures.tsv` records each workload's role, cluster and protein counts, and
archive checksum. The extraction process checks the checksum before unpacking.

## Versions and image

The Dockerfile pins Python Clinker 0.0.32 and Rust `clinker-cli` 0.1.0, plus
the Python dependency versions in `requirements.txt`. It builds the Rust CLI
with Rust 1.88 in a separate stage, then places the commands at
`/opt/python-clinker/bin/clinker` and `/usr/local/bin/clinker-rust`. Rust 1.87
is needed by the CLI's locked `bio` dependency, and Rust 1.88 supports the
edition-2024 let chains and `slice::as_chunks` used by the published CLI. The
base images are pinned by manifest-list digest. Build and publish the image
once, then use its immutable digest for release measurements. For local
rebuilds, use a new image tag when using `-resume`; otherwise Nextflow could
reuse completed tasks from the previous image with the same tag.

The workflow caps benchmark concurrency at one task to avoid CPU contention
between measured runs. Nextflow can still resume successful tasks after a
failure as long as its cache and work directory remain available.
