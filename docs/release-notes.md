# Initial Rust release

- Browser-based GenBank and GFF3/FASTA comparison with local WASM analysis.
- Cancellable analysis, optional k-mer prefilter, bounded Ward ordering, and
  automatic locus orientation and placement.
- Plot editing, SVG and similarity CSV export, and complete saved projects.
- Native CLI builds for Linux x86-64, Windows x86-64 and Apple Silicon macOS.
- Embedded local viewer with tree/annotation import, editing and exports; native
  results are rendered without WASM analysis. Bare `--plot` opens the viewer.
- Shared editable project format; `--project` writes a complete project and
  `--open-project` restores it with its figure edits.
- Self-contained core/CLI Cargo source packages with independent install checks,
  native archive extraction checks and checksums for every release artifact.

The Rust port does not reproduce every Python score or layout exactly. Python
session compatibility and adding input files to
an edited browser project are outside this release. Similarity CSV scores come
from the retained analysis; visual edits do not rerun it. Translation fallback
does not infer alternative initiation or recoding exceptions.

Approximate prefiltering is opt-in in both the browser and CLI. Its performance
benefit depends on the input, and it can discard distant homologues.

CLI and browser projects retain the computed Ward dendrogram when similarity
ordering runs. Choose Row order → Similarity dendrogram to display it; saved
tree edits can be restored in either viewer. See the [acceptance baseline and
remaining release gates](features.md#local-acceptance-baseline--2026-10-27).
