# Initial Rust release

- Browser-based GenBank and GFF3/FASTA comparison with local WASM analysis.
- Cancellable analysis, optional k-mer prefilter, bounded Ward ordering, and
  automatic locus orientation and placement.
- Plot editing, SVG and similarity CSV export, and complete saved projects.
- Native CLI builds for Linux x86-64, Windows x86-64 and Apple Silicon macOS.

The Rust port does not reproduce every Python score or layout exactly. Python
session compatibility, a CLI-hosted interactive app, and adding input files to
an edited browser project are outside this release. Similarity CSV scores come
from the retained analysis; visual edits do not rerun it. Translation fallback
does not infer alternative initiation or recoding exceptions.
