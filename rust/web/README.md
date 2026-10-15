# Static web prototype

This is a deliberately small static frontend. It accepts GenBank files and
GFF3 annotations with a same-stem companion FASTA or embedded `##FASTA` section. The browser sends uploaded files
to `worker.js` once for parsing. It then sends small protein-block tiles to a
pool of `worker.js` instances; the workers perform up to 1,600 gene comparisons
per tile and return retained links. Those links are combined with the original
cluster metadata and drawn by the existing `clustermap.js` assets. No file
contents are uploaded.

The plot's data and appearance controls use clustermap.js's optional
`<clinker-editor>` Web Component. The build copies its ESM entry point and
default stylesheet alongside the main renderer, so the editor remains fully
offline-capable on static hosting.

During a tile, WASM reports progress every 100 processed pairs. The page sums
those counters across workers, so the analysis breadcrumb reports protein
pairs rather than completed tiles. Parsing, candidate filtering, cluster-pair
scoring, homology grouping, and locus arrangement also expose their own
measured work units through the shared Rust core.

The optional **k-mer prefilter** indexes distinct protein 3-mers and only
aligns pairs with at least three shared words. Candidate tiles contain each
needed protein once plus pair indices, rather than one copied sequence per
pair. This is deliberately approximate and can miss remote homologues. It is
enabled by default with permissive settings (k-mer size 3, three shared words);
disable the checkbox beside the identity setting to compare every eligible
protein pair.
Larger words and higher shared-word thresholds are more selective. Protein
normalisation is shared with alignment. Words are cached once per protein and
only one target-cluster index is held at a time. Words use lossless five-bit
residue encoding: 32-bit keys for k ≤ 6 and 64-bit keys for k = 7–8. Cached word
lists discard unused capacity after deduplication. Temporary caches are released
with the parsing worker before alignment begins.
A separate filtering step reports the number and percentage
of possible alignments removed before the remaining candidates are aligned.
Unfiltered mode still skips protein pairs whose length ratio
makes the chosen global-identity cutoff mathematically impossible.
The run panel exposes the k-mer size and shared-word threshold as advanced
settings, along with the choice between synteny ordering and input file order.
It can also select one-based inclusive record ranges, split multi-record files
into separate clusters, and opt out of automatic circular-origin rotation.
Explicit trim intervals always refer to the original record coordinates and
suppress automatic origin rotation for that record.
Default row ordering builds a Ward hierarchy, then approximately minimizes
neighbouring distances while preserving every hierarchical group. Leaf ordering
retains at most 32 endpoint alternatives per subtree, storing reconstruction
pointers rather than full candidate sequences. With this fixed search width,
leaf ordering takes linear time and memory; the Ward hierarchy and pairwise
matrices use O(n²) time and memory. The ordering is deterministic but does not
guarantee the optimum among all tree-compatible orders.
Annotation-only runs can skip protein alignment entirely. Optional two-column
gene-function and function-colour CSV files seed labelled, coloured groups;
retained homology links extend those groups when alignment is enabled.

Repository test examples are available in
[`gene-functions.csv`](../../examples/gene-functions.csv) and
[`function-colours.csv`](../../examples/function-colours.csv) files. They seed
three illustrative groups using protein IDs in the bundled *A. versicolor*
cluster and work with either its GenBank or GFF3 input, including **Try example**.
These are demonstration labels, not biological function annotations. CSV files
have no header row. With alignment enabled, each seeded group extends to linked
genes in the other example clusters.

The web app has one restorable export: a versioned clinker project. It contains
the plot data, renderer appearance, editing state, and complete analysis state,
including translations, links, identity cutoff, and configured function groups.
This versioned project is the only file accepted by the web app's Open action.
Reopening restores the saved figure directly and computes similarity scores
from its retained analysis; it does not rerun ordering, grouping, or layout.
Similarity CSV exports use current cluster names and plot positions. Scores
remain those of the original analysis, including after genes are edited or
deleted from the figure.
Adding new files to an existing analysis is currently available through the
CLI; browser append support remains pending.

GFF3 annotations and companion FASTA files stay together in the input-file
list and are removed together. Counts include every selected input file.
The prefilter checkbox sits beside the identity setting; its k-mer size and
shared-word threshold are in Advanced settings.
Analysis, example downloads, and project loading can be cancelled. Cancelling
terminates their workers and keeps any existing plot and selected inputs.
The progress view retains the analysis stages and their progress bars.
Input pairing and trim controls share one sequential inspection pass. Only
sequence metadata is cached; changing the selection aborts obsolete reads.

Run `npm test` for input pairing, workflow lifecycle, and project/export tests.
After building WASM, run `npm run test:wasm` to exercise the incremental parser
and exported analysis functions. For the browser round trip, install Chromium
with `npx playwright install chromium` and run `npm run test:browser`. CI runs
that browser test against the packaged static application. It also checks
cancellation at worker-stage boundaries, worker failures, and retrying project
loading while preserving an edited plot. The editor round trip exercises bulk
gene deletion through a locus header, undo/redo, and saving the resulting figure.

## Build for local preview or static hosting

The development machine needs Node/npm, Rust's `wasm32-unknown-unknown` target,
and [`wasm-pack`](https://rustwasm.github.io/docs/wasm-pack/). `npm ci` fetches
the pinned clustermap.js refactor commit and D3 v7 from `package-lock.json`.
The Homebrew Rust installation can coexist with rustup, but the build script
explicitly selects rustup's Cargo toolchain because that is where the WASM
target is installed.

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack
cd rust/web
sh build.sh
python3 -m http.server
```

`build.sh` finds Cargo-installed `wasm-pack` at `~/.cargo/bin/wasm-pack`, so
adding that directory to your shell `PATH` is optional.

Release builds currently skip wasm-opt because wasm-pack does not provide a
usable Apple Silicon binary in this setup. This changes file size, not analysis
correctness.

Open the reported local URL; do not open `index.html` directly because module
workers and WASM modules require HTTP serving. Deploy the contents of this
directory, including `pkg/`, to GitHub Pages or another static host.
`sh package.sh /path/to/site` copies only runtime assets and examples, excluding
test dependencies and sources. Pages and release archives use this same package.
