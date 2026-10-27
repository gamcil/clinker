# Feature coverage and release scope

The native CLI and browser share the Rust analysis implementation. Supported
workflows are listed below; scores, traceback choices and cluster order are not
required to match the legacy Python implementation exactly.

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
| Display the computed Ward dendrogram | Shared renderer tree | Local viewer and saved projects | Yes | Core hierarchy, CLI project and browser tree regressions |
| Offline analysis without Python or an external service | Yes | Yes | Yes | Native and WASM builds |

Newick and annotation imports belong to the JavaScript host. Rust project
persistence preserves their renderer data as JSON.

## Replacement workflow status

The native `--plot` option now opens the embedded viewer when no path is given.
Packaging is implemented. Remaining release work includes cross-platform and
hosted verification and the deferred repository cutover.

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

## Data imports

Use **Add data** to import Newick trees or TSV/CSV annotations. Tree leaves must
match every plotted cluster exactly once; ambiguous names can be mapped manually.
Quoted Newick labels, non-negative branch lengths and internal support values
are supported. Trees can be rotated and styled, and edits survive project saving.
**Row order** switches between the default order, imported tree and computed
similarity dendrogram. The latter represents gene-cluster similarity, not ancestry.

For annotations, **Separate tracks** gives each selected column its own scale
or legend and uses its header as the panel title. Use **One heatmap** for columns
with comparable units, such as expression values or presence/absence. It shares
one scale; values are not automatically normalized. Blank cells remain missing.
See the [example table](../examples/annotations.md).

## Release checks

CI checks native Linux/macOS/Windows builds, source packages and the packaged web
app. Manual verification of downloaded releases on supported systems and the
hosted app in actual Safari remains part of the [release procedure](releasing.md).
