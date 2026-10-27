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
| Seed groups from a gene-function CSV | Yes | Yes | Yes | Group and CLI tests; WASM build |
| Apply function colours from CSV | Yes | Yes | Yes | Group and CLI tests; WASM build |
| Synteny-based ordering or input order | Yes | Yes | Yes | Ordering tests |
| Automatic locus arrangement and orientation | Yes | Via plots | Yes | Layout and plot-data tests |
| Human-readable or delimited link summary | Structured links | Yes | Not applicable | CLI workflow |
| Similarity/distance matrix export | Yes | CSV | CSV | Matrix tests and CLI smoke test |
| clustermap.js plot-data output | Yes | JSON | In project | Plot-data tests |
| Portable interactive figure | Yes | HTML and local viewer | Static application | CLI and browser smoke tests |
| Save and reopen complete analysis state | Yes | Session; project export | Project | Core session round-trip, CLI append workflow, and browser smoke test; native project/viewer flow checked separately |
| Add new files to a saved session | Yes | Yes | Future enhancement | Core and CLI workflow tests |
| SVG export and interactive plot editing | Plot contract | Local viewer | Yes | Browser integration coverage and native viewer checks |
| Import Newick trees and preserve interactive tree edits | Tree data | Local viewer | Yes | Parser/matching fixtures and browser project/export tests |
| Display the computed Ward dendrogram | Hierarchy | Not yet attached to newly generated native projects; restored from saved projects | Yes | Browser tree tests |
| Offline analysis without Python or an external service | Yes | Yes | Yes | Native and WASM builds |

## Replacement workflow status

The native `--plot` option now opens the embedded viewer when no path is given.
Remaining migration work concerns distribution packaging and release cutover.

## Future additions

| Enhancement | Current state | Release scope |
| --- | --- | --- |
| Add new data to an edited browser project | Core and CLI can append to saved analysis sessions, and the web app can reopen and export projects. A browser workflow would also need to reconcile new analysis results with existing figure edits. | Post-1.0; not required for Python parity |
| Import files directly from external web apps | Deferred until the web app is more mature and there is concrete integration demand. A future handoff would populate input files in a new clinker tab, leaving settings and starting analysis to the user. | Later addition; not required for release or Python parity |

## Release criterion

Each claimed feature should have a fixture-backed core test and, where exposed,
a CLI or browser workflow test. A release does not require output equivalence
with Python, but unsupported established workflows must remain visible in this
matrix and in release notes.
