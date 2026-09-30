# Static web prototype

This is a deliberately small static frontend. The browser sends uploaded files
to `worker.js` once for parsing. It then sends small protein-block tiles to a
pool of `worker.js` instances; the workers perform up to 1,600 gene comparisons
per tile and return retained links. Those links are combined with the original
cluster metadata and drawn by the existing `clustermap.js` assets. No file
contents are uploaded.

During a tile, WASM reports progress every 100 processed pairs. The page sums
those counters across workers, so the status line reports protein pairs rather
than completed tiles.

The optional **Fast k-mer prefilter** indexes distinct protein 3-mers and only
aligns pairs with at least three shared words. Candidate tiles contain each
needed protein once plus pair indices, rather than one copied sequence per
pair. This is deliberately approximate and can miss remote homologues, so it
is off by default. Exact mode still skips protein pairs whose length ratio
makes the chosen global-identity cutoff mathematically impossible.

## Build for local preview or static hosting

The development machine needs Rust's `wasm32-unknown-unknown` target and
[`wasm-pack`](https://rustwasm.github.io/docs/wasm-pack/). The Homebrew Rust
installation can coexist with rustup, but the build script explicitly selects
rustup's Cargo toolchain because that is where the WASM target is installed.

```bash
rustup target add wasm32-unknown-unknown
cargo install wasm-pack
cd rust/web
sh build.sh
python3 -m http.server
```

`build.sh` finds Cargo-installed `wasm-pack` at `~/.cargo/bin/wasm-pack`, so
adding that directory to your shell `PATH` is optional.

The WASM package explicitly enables `getrandom`'s browser backend because the
shared `bio` dependency reaches it transitively. No browser-side randomness is
used by clinker itself.

Release builds currently skip wasm-opt because wasm-pack does not provide a
usable Apple Silicon binary in this setup. This changes file size, not analysis
correctness.

Open the reported local URL; do not open `index.html` directly because module
workers and WASM modules require HTTP serving. Deploy the contents of this
directory, including `pkg/`, to GitHub Pages or another static host.
