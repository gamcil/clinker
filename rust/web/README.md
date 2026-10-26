# Static web prototype

This is a deliberately small static frontend. It accepts GenBank files and
GFF3 annotations with a same-stem companion FASTA or embedded `##FASTA` section. The browser sends uploaded files
to `worker.js` once for parsing. It then sends small protein-block tiles to a
pool of `worker.js` instances; the workers perform up to 1,600 gene comparisons
per tile and return retained links. Those links are combined with the original
cluster metadata and drawn by the existing `clustermap.js` assets. No file
contents are uploaded.

The plot's data and appearance controls use clustermap.js's optional
`<clinker-editor>` Web Component. The build copies its ESM entry point, optional
annotation overlay and shared module chunks alongside the main renderer. Styles
are embedded in those modules, so the editor remains fully offline-capable on
static hosting. No registration call or separate editor stylesheet is needed.

## Application and editor styling

Clinker owns its application controls, dialogs and import workflows. App markup
is scoped under `.clinker-app`; ordinary controls use `.app-button`, `.app-field`
and `.app-disclosure`, including controls created dynamically. `style.css` defines
semantic `--app-*` values for both palettes and shares component rules between
light and dark modes. Primary/status actions have explicit foreground/background
pairs. Keep new app styles inside this scope and use these variables instead of
adding another block of dark-mode component overrides.

Clustermap owns editor internals, selection actions and floating plot controls.
The `.clinker-plot-workspace` adapter maps app tokens to its public `--cm-*`
properties for the chart, editor and selection toolbar. Keep library-specific
rules limited to host placement and dimensions; avoid selectors for its internal
buttons, borders or text. Do not query editor shadow roots for navigation.

The dependency is pinned to refactor commit
`f532521cebf919222d9b7865b54a1c3c5af6b0d7`. Clinker sets the application's saved
System/Light/Dark preference with `chart.uiAppearance(...)`; loading a project
or editing figure settings does not change that preference. Interface appearance
does not modify figure colours, annotation scales, undo history or exported SVG.
Tree appearance uses `editor.openAppearance({ section: "tree", focus: true })`.
The public `clinker-editor-open` and `clinker-editor-close` events synchronise
the toolbar button, including when selection editing opens the editor.

Build/package scripts copy the core bundle, `editor.mjs` (renamed locally to
`clustermap-editor.mjs`), `annotationOverlay.mjs`, and the entire `editor/` chunk
directory. Preserve their relative module paths. The CLI's vendored core bundle
must be updated alongside the browser dependency.

Independent annotation tracks use the panel title as their single heading and
hide redundant column labels. A multi-column heatmap retains column labels and
can have an optional overall title. Automatic spacing of wide legends remains
an upstream issue documented in `CLUSTERMAP_ANNOTATION_LAYOUT_HANDOFF.md` at the
repository root; the UI-boundary update does not correct that layout algorithm.

## Analysis workflow

During a tile, WASM reports progress every 100 processed pairs. The page sums
those counters across workers, so the analysis breadcrumb reports protein
pairs rather than completed tiles. Parsing, candidate filtering, cluster-pair
scoring, homology grouping, and locus arrangement also expose their own
measured work units through the shared Rust core.

Supplied protein translations take precedence. Translation fallback honours
GenBank `codon_start`, GFF CDS phase, and `transl_table` in either format using
the [NCBI sense-codon tables](https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi).
Initial frame offsets apply after splicing, on the biological strand. Alternative
initiation, `transl_except`, and context-dependent recoding are not inferred;
use an annotated protein translation for those cases.

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
Locus orientation scores collinear homology links against linked, already placed
loci. Scoring reuses one DP row; placement uses weighted-median anchors without
traceback. Unanchored loci retain deterministic fallback positions.

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

## Trees

Choose an optional **Tree** directly in the new-run panel, or use
**Order → Import tree** in the plot toolbar after analysis. New-run tree matches
appear below the file chooser and edits take effect immediately; no separate
confirmation is needed. Imports into an existing plot use a dialog. A tree selected for a new run is matched using the inspected input filenames
and sequence IDs, so differences can be resolved before starting. Matches are
checked against parsed clusters before protein alignment; if the parsed inputs
differ, analysis pauses for corrections. The computed default order is still
retained. The
dialog matches each leaf to exactly one plotted cluster and lists any missing,
ambiguous or duplicate matches. It accepts exact cluster names, cluster UIDs,
locus/sequence names, source filenames, basenames and annotation filename
stems. Matching is case-sensitive: accession versions, underscores and spaces
are preserved. Use the mapping dropdowns to resolve differences without
renaming the inputs. A multi-locus cluster is one tree leaf; records plotted
as separate clusters each require their own leaf. Extra leaves and missing
clusters must be resolved before applying; automatic pruning is not supported.

The reader accepts one semicolon-terminated Newick tree, single-quoted labels
with doubled apostrophes (`'O''Brien'`), optional finite non-negative branch
lengths (including zero and scientific notation), named internal nodes,
multifurcations and unary nodes. Quote labels containing spaces or punctuation.
Leading `[&R]` and `[&U]` flags set biological rootedness; otherwise it remains
unspecified. Other comments and annotations, including NHX, are rejected with
an error rather than discarded. The parser is a small iterative frontend
module; no tree-viewer or parser dependency is added.

Under **Import options**, choose whether numeric internal labels are support
values or node names. Numeric leaf labels remain leaf labels. Original internal text is kept
even when also interpreted as support. Support scale defaults to unspecified;
percent and proportion are explicit alternatives and values are never rescaled.
Incoming branch lengths/support on the display root require the **Discard root branch length and support** option before import, shown only when needed. This discards those root attributes
while retaining the root label. Missing branch lengths are retained as unknown;
the viewer falls back to a cladogram if any non-root length is absent or every
root-to-tip distance is zero.

**Order → Show clustering dendrogram** displays the Ward hierarchy already used
for similarity ordering. Its initial leaf order matches the computed plot order;
branches can be rotated using the same controls as an imported tree. It shows
gene-cluster similarity, with horizontal spacing scaled by Ward merge height
and all leaves aligned at zero height. **Scale by merge height** in the Order
menu toggles between this default and an unscaled topology. The same setting
is available as **Use branch lengths** under tree appearance. These are
clustering distances, with no support values.
The dendrogram and any imported tree are retained separately when switching
orders, and both survive project saving and reopening. It is available for new
analyses that perform similarity ordering; keeping input order or skipping
alignment does not generate a dendrogram. Older projects without a saved
hierarchy continue to open normally.

The [bundled example tree](../../examples/bundled-clusters.nwk) matches **Try
example** and is fully bifurcating. Its topology, branch lengths and support values are illustrative,
not an inferred phylogeny.

Click an internal node to reverse its children. Vertical cluster and batch
dragging follows the tree's clade constraints; horizontal locus movement,
flipping and trimming remain independent. Node rotations and compatible drags
support undo/redo. Use the existing **Edit plot → Appearance → Tree** section to
change width, branches, support labels, guides and the tree scale bar.
Right-click an internal node for branch/clade styling and text annotations;
these edits support undo/redo and
are preserved in saved projects and SVG exports. Generated Ward dendrograms
also expose cut-height clade colouring in the Tree appearance controls. Hiding
the tree keeps its ordering constraints. The toolbar's **Order** menu
switches between **Computed** and **Tree**. The computed default is the
original analysis order retained in the plot-data cluster array; tree rotations
and manual dragging change display state, not that array. The menu shows
**Manual** after unrestricted manual reordering, so selecting Computed
also resets those edits. Switching to computed order retains the edited tree
for switching back. **Order → Remove tree** discards it instead.
Applying/replacing/removing a tree or switching tree ordering resets chart undo
history; existing gene edits and locus placements are retained.

Version 2 projects preserve the live edited tree, child order, appearance and
analysis data. When using computed order, the retained tree is saved in the
optional frontend `inactiveTree` field and validated on reopening. The Ward
dendrogram is saved separately as `clusteringDendrogram`, with `treeSource`
identifying when it is active. Invalid tree projects are rejected before
replacing the figure.
Tree file reads can be cancelled by closing the dialog; new analyses/project
loads cancel pending imports. Each new analysis rematches the tree selected
in the new-run panel, rather than reusing the previous figure's index-based
mapping. Clear that selection to run without a tree. SVG export includes the tree in every
renderer. Tree inference and rerooting are not provided.

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
with `npx playwright install chromium` and run `npm run test:browser`. Set
`CLINKER_TEST_BROWSER=firefox` or `webkit` to use another installed Playwright
browser. CI installs all three engines and runs the test against the packaged
static application in each. It also checks
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
