# Feature coverage and release scope

The Rust implementation is intended to become clinker, rather than remain a
second implementation that reproduces every numerical result from the Python
package. The Python application is therefore used to identify established user
workflows. Rust tests verify those capabilities directly with fixtures and
end-to-end workflows.

Alignment scores, traceback choices, link counts near a cutoff, cluster order,
and generated identifiers do not need to match Python exactly. They must be
well-defined, useful, stable enough for saved output, and covered by Rust tests.

## Implemented

| Capability | Core | CLI | Web | Verification |
| --- | --- | --- | --- | --- |
| Parse one or more GenBank files and multi-record files | Yes | Yes | Yes | Parser fixtures and bundled workflow |
| Parse GFF3 with companion or embedded FASTA | Yes | Yes | Yes | Core fixtures, CLI workflow, and WASM adapter test |
| Discover supported annotation files recursively | Not applicable | Yes | Not applicable | CLI discovery test |
| Select `record:start-end` ranges | Yes | Yes | Yes | Parser, CLI, and WASM adapter tests |
| Treat records in one file as separate clusters | Yes | Yes | Yes | Parser and WASM adapter tests |
| Retain source-file and GenBank organism provenance | Yes | Via plot data | Yes | Bundled and separate-record parser fixtures |
| Rotate features crossing a circular origin, with an opt-out | Yes | Yes | Yes | Circular GenBank fixture |
| Forward, reverse, translated, derived, and compound CDS features | Yes | Yes | Yes | Parser fixtures |
| Cross-cluster global protein alignment with an identity cutoff | Yes | Yes | Yes | Alignment and workflow tests |
| Control native alignment concurrency | Local Rayon pool | `-j/--jobs` | Web Workers | Serial/parallel equality and CLI session tests |
| Parse and plot without protein alignment | Yes | Yes | Yes | Annotation-only group tests |
| Single-linkage homology groups | Yes | Via plot data | Yes | Group and plot-data tests |
| Seed groups from a gene-function CSV | Yes | Yes | Yes | Group/CLI regressions and WASM group round-trip |
| Apply function colours from CSV | Yes | Yes | Yes | Group/CLI regressions and WASM group round-trip |
| Synteny-based ordering or input order | Yes | Yes | Yes | Ordering tests |
| Automatic locus arrangement and orientation | Yes | Via plots | Yes | Layout and plot-data tests |
| Human-readable or delimited link summary | Structured links | Yes | Not applicable | CLI workflow |
| Similarity/distance matrix export | Yes | CSV | CSV | Matrix tests and CLI smoke test |
| clustermap.js plot-data output | Yes | JSON | In project | Plot-data tests |
| Portable interactive figure | Yes | HTML and local viewer | Static application | CLI and browser smoke tests |
| Save and reopen complete analysis state | Yes | Session; project export | Project | Core session round-trip, CLI append workflow, and browser smoke test; native project/viewer flow checked separately |
| Add new files to a saved session | Yes | Yes | Future enhancement | Core and CLI workflow tests |
| SVG export and interactive plot editing | Plot contract | Local viewer | Yes | Browser integration coverage and native viewer checks |
| Import TSV/CSV annotation tracks and shared-scale heatmaps | Project round-trip | Local viewer | Yes | Browser import/edit/project round-trip; native viewer acceptance |
| Import Newick trees and preserve interactive tree edits | Project round-trip | Local viewer | Yes | Parser/matching fixtures and browser project/export tests |
| Display the computed Ward dendrogram | Hierarchy | Not yet attached to newly generated native projects; restored from saved projects | Yes | Browser tree tests |
| Offline analysis without Python or an external service | Yes | Yes | Yes | Native and WASM builds |

Newick and annotation imports belong to the JavaScript host. Rust project
persistence preserves their renderer data as JSON.

## Replacement workflow status

The native `--plot` option now opens the embedded viewer when no path is given.
Packaging is implemented. Remaining release work includes cross-platform and
hosted verification, the deferred repository cutover, and deciding whether the
native computed-dendrogram gap above must be closed before the first release.

## Future additions

| Enhancement | Current state | Release scope |
| --- | --- | --- |
| Add new data to an edited browser project | Core and CLI can append to saved analysis sessions, and the web app can reopen and export projects. A browser workflow would also need to reconcile new analysis results with existing figure edits. | Post-1.0; not required for Python parity |
| Import files directly from external web apps | Deferred until the web app is more mature and there is concrete integration demand. A future handoff would populate input files in a new clinker tab, leaving settings and starting analysis to the user. | Later addition; not required for release or Python parity |

## Release criterion

Each claimed workflow needs a substantive regression at its owning layer and a
representative CLI or browser integration check where exposed. Follow the
[test ownership policy](testing.md); detailed renderer/editor checks belong to
clustermap. A release does not require output equivalence with Python, but
unsupported established workflows must remain visible in this matrix and in
release notes.

## Local acceptance baseline — 2026-10-27

Measured on macOS 26.5.2, arm64, eight logical CPUs, using Rust 1.98.1 for both
native and WASM, Node 25.8.2 and Chromium 153.0.8010.12. Each configuration ran
three times with four workers and identity cutoff 0.30. The optional filter used
k=3 and three shared words. Native measurements used the release-mode core
benchmark; browser measurements used the existing Chromium/Web Worker harness.
Input staging and module compilation are excluded. Timings include parsing,
preparation/filtering and alignment, and exclude ordering, layout and rendering.

The small input is the five frozen GenBank examples: 51 proteins and 1,001
possible comparisons. The larger local `bgcs_gcf_export_2` snapshot has 169
clusters, 1,465 proteins and 1,065,430 possible comparisons. Input snapshot
fingerprints are `2863f5a436df` and `e97ad6ab38ec`, respectively (SHA-256 over
sorted filename/NUL/file-SHA-256/newline records). The larger dataset is not
part of the distributed fixtures.

| Input | Filter | Native median (s) | WASM median (s) | Retained links, both |
| --- | --- | --- | --- | --- |
| Five examples | Off | 0.034 | 0.108 | 79 |
| Five examples | 3/3 | 0.095 | 0.170 | 79 |
| 169-cluster GCF export | Off | 22.263 | 22.714 | 6,784 |
| 169-cluster GCF export | 3/3 | 23.602 | 24.839 | 6,778 |

Native/WASM protein, candidate, alignment and retained-link counts agreed on every
run. On the larger input the unfiltered path ran 648,549 alignments after the
exact length bound; the filter selected 611,618 pairs and retained six fewer
links. Filtering was slower on both measured inputs. Filtering remains opt-in in
both frontends; the browser previously enabled it by default.

Native peak RSS on the larger input was 34.0 MiB unfiltered and 45.9 MiB filtered.
Sampled Node/Chromium process-tree RSS reached approximately 960 and 913 MiB.
Those browser totals include browser overhead and double-count shared pages;
they are not unique physical-memory measurements or directly comparable to the
native process. Larger inputs remain unmeasured.

A full unfiltered CLI run on the larger input, writing a session, editable
project, portable HTML, link summary and distance CSV, took 22.84 seconds
and peaked at 71.6 MiB RSS. The packaged browser app took 27.31 seconds
from clicking Analyse to displaying the same 169-cluster, 6,784-link figure
(single end-to-end observations).

Verified locally:

- Packaged application workflows in Chromium, Firefox and WebKit, including
  separate TSV tracks, a multi-column CSV presence/absence heatmap and project
  persistence through edits and reopening.
- The extracted optimized macOS archive: native viewer, editing, annotation
  import, project reopening, similarity CSV and SVG export on the larger input.
- Portable HTML rendering directly from disk.
- Extracted core/CLI Cargo archives: build, installation, exports and every
  embedded viewer asset under Cargo 1.98. The package checker now initializes
  its dedicated cache marker before cleaning.

Remaining release gates:

- Fresh CLI projects retain computed row order but omit the Ward dendrogram.
  Saved browser-generated dendrograms are restored correctly. Resolve this
  entry-point difference or explicitly accept it for the first release.
- Confirm Linux/Windows native archives, actual Safari and the hosted Pages
  deployment in their release environments, following [releasing.md](releasing.md).
  Local WebKit checks cover engine behavior; Safari verification remains pending.
