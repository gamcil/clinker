# Embedded viewer assets

`viewer/` is a complete, generated snapshot of the shared local viewer: HTML,
application controls, styles, renderer, editor modules, hashed chunks and licence
notices. Cargo embeds this committed snapshot using the Rust-only `build.rs`;
`cargo build` and `cargo install` never run npm, wasm-pack or a frontend download.

The native viewer starts no WASM analysis workers. It reads projects and
similarity information from the CLI's loopback server. The portable HTML exporter
also embeds its D3 and renderer bundles from this snapshot, without the editor.

Source versions are pinned by `web/package-lock.json`. From the repository root:

```sh
npm --prefix web ci --ignore-scripts
npm --prefix web run sync:assets
```

This builds the viewer from shared source without compiling WASM, replaces the
snapshot completely, and synchronizes the root licence notices. Commit all
changes, including renamed or removed content-hashed chunks. The source inventory
is `web/scripts/assets.mjs`; `web/scripts/build-viewer.mjs` selects the shared
shell and local entry. The viewer's build ID excludes WASM, uploaded/example
inputs and the browser analysis entry.

CI runs `npm --prefix web run check:assets`. It rebuilds the expected viewer and
compares both file inventory and bytes without changing committed assets.
Normal web builds do not rewrite this snapshot.
