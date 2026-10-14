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
| clustermap.js plot-data export | Yes | JSON | JSON | Plot-data tests |
| Portable interactive figure | Yes | HTML | Static application | CLI and browser smoke tests |
| SVG export and interactive plot editing | Plot contract | Not applicable | Yes | Browser integration coverage |
| Offline analysis without Python or a server | Yes | Yes | Yes | Native and WASM builds |

## Features still needed for replacement

| Established workflow | Current gap | Suggested priority |
| --- | --- | --- |
| Save and resume an analysis session | No stable analysis-session schema | Medium |
| Add new files to a saved session | Depends on the session model | Medium |
| Control native alignment concurrency | Native core is serial; web uses workers | Medium |
| Open an interactive plot when `--plot` has no path | CLI currently requires an output path | Low; the static web app covers this workflow |

## Release criterion

Each claimed feature should have a fixture-backed core test and, where exposed,
a CLI or browser workflow test. A release does not require output equivalence
with Python, but unsupported established workflows must remain visible in this
matrix and in release notes.
