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
| Parse and plot without protein alignment | Yes | Yes | Yes | Annotation-only group tests |
| Single-linkage homology groups | Yes | Via plot data | Yes | Group and plot-data tests |
| Seed groups from a gene-function CSV | Yes | Yes | Yes | Group and CLI tests; WASM build |
| Apply function colours from CSV | Yes | Yes | Yes | Group and CLI tests; WASM build |
| Synteny-based ordering or input order | Yes | Yes | Yes | Ordering tests |
| Automatic locus arrangement and orientation | Yes | Via plots | Yes | Layout and plot-data tests |
| Human-readable or delimited link summary | Yes | Yes | Not applicable | Formatter and CLI tests |
| Similarity/distance matrix export | Yes | CSV | CSV | Matrix tests and CLI smoke test |
| clustermap.js plot-data output | Yes | JSON | In project | Plot-data tests |
| Portable interactive figure | Yes | HTML | Static application | CLI and browser smoke tests |
| Save and reopen complete analysis state | Yes | Session | Project | Core round-trip, CLI append workflow, and browser smoke test |
| Add new files to a saved session | Yes | Yes | Future enhancement | Core and CLI workflow tests |
| SVG export and interactive plot editing | Plot contract | Not applicable | Yes | Browser integration coverage |
| Offline analysis without Python or a server | Yes | Yes | Yes | Native and WASM builds |

## Features still needed for replacement

| Established workflow | Current gap | Suggested priority |
| --- | --- | --- |
| Control native alignment concurrency | Native core is serial; web uses workers | Medium |
| Open an interactive plot when `--plot` has no path | CLI currently requires an output path | Low; the static web app covers this workflow |

## Future additions

| Enhancement | Current state | Release scope |
| --- | --- | --- |
| Add new data to an edited browser project | Core and CLI can append to saved analysis sessions, and the web app can reopen and export projects. A browser workflow would also need to reconcile new analysis results with existing figure edits. | Post-1.0; not required for Python parity |

## Release criterion

Each claimed feature should have a fixture-backed core test and, where exposed,
a CLI or browser workflow test. A release does not require output equivalence
with Python, but unsupported established workflows must remain visible in this
matrix and in release notes.
