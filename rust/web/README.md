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
pair. This is deliberately approximate and can miss remote homologues, so it
is off by default. A separate filtering step reports the number and percentage
of possible alignments removed before the remaining candidates are aligned.
Unfiltered mode still skips protein pairs whose length ratio
makes the chosen global-identity cutoff mathematically impossible.
The run panel exposes the k-mer size and shared-word threshold as advanced
settings, along with the choice between synteny ordering and input file order.
It can also select one-based inclusive record ranges, split multi-record files
into separate clusters, and opt out of automatic circular-origin rotation.
Annotation-only runs can skip protein alignment entirely. Optional two-column
gene-function and function-colour CSV files seed labelled, coloured groups;
retained homology links extend those groups when alignment is enabled.

The web app has one restorable export: a versioned clinker project. It contains
the plot data, renderer appearance, editing state, and complete analysis state,
including translations, links, identity cutoff, and configured function groups.
This versioned project is the only file accepted by the web app's Open action.
Adding new files to an existing analysis is currently available through the
CLI; browser append support remains pending.

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
