# Architecture

The Cargo workspace contains three packages sharing one lockfile:

| Package | Responsibility |
| --- | --- |
| `clinker-core` | Parsing bytes, alignment, filtering, grouping, ordering, layout and saved analysis |
| `clinker-cli` | File discovery, native execution, exports, persistence and the local viewer |
| `clinker-wasm` | JavaScript conversion and worker-facing bindings to the core |

Core records live in `model.rs`; alignment options, batches and worker results
live in `alignment.rs`. Native Rayon threads and browser Web Workers execute
the same batch planner and alignment loop, reusing buffers. Unfiltered pairs
are expanded within batches; optional k-mer filtering stores candidate pairs.

`Analysis::prepare_figure` scores cluster pairs, orders rows, groups genes and
arranges loci. It returns plot data, input-indexed similarities, ordering and a
renderer-ready Ward dendrogram. `RowOrder::Input` skips similarity ordering and
the dendrogram. `to_plot_data(order, definitions)` only converts existing
analysis, preserving locus coordinates and placement.

Analysis retains translations and original records. Renderer data uses stable
IDs and omits translations. `AnalysisSession` supports native continuation;
`ClinkerProject` also retains figure edits. See [project formats](project-format.md).

## Browser and native viewer

`web/src/app.js` owns browser inputs and analysis orchestration; `local.js`
loads native results. Both use `viewer/` for chart lifecycle, menus, appearance
and exports. Tree and annotation parsing/matching belong to Clinker; rendering,
editor controls and figure interactions belong to clustermap.

The native shell is built from shared HTML with run controls omitted.
`crates/clinker-cli/assets/viewer/` is the committed generated snapshot, embedded
by a Rust-only build script. Native consumers need no Node.js or WASM. Portable
HTML is a separate renderer-only export. See [packaging](packaging.md).

## UI integration

App controls use `.app-button`, `.app-field` and `.app-disclosure` under
`.clinker-app`, with semantic `--app-*` theme variables. The
`.clinker-plot-workspace` adapter maps these to clustermap's public `--cm-*`
properties. Keep selectors out of library internals and editor shadow roots.

Use `chart.uiAppearance(...)` for System/Light/Dark preferences, independently
of figure colours or saved projects. Use
`editor.openAppearance({ section: "tree", focus: true })` to navigate to tree
settings. Public `clinker-editor-open` / `clinker-editor-close` events keep host
controls synchronized. Preserve upstream module paths and hashed chunk names.

## Development

```sh
cargo test --locked --workspace
cargo clippy --locked --workspace --all-targets -- -D warnings
cargo fmt --all -- --check
```

See [web development](../web/README.md), [test ownership](testing.md),
[CLI usage](cli.md) and [release procedure](releasing.md).
