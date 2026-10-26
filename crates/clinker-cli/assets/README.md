# Embedded plot assets

The CLI embeds `d3.min.js` and `clustermap.min.js` to make `--plot *.html`
self-contained. These generated bundles are deliberately vendored with the
crate so native Cargo builds do not depend on Node, npm, or network access.

Their source versions are pinned by `web/package-lock.json`. From the repository
root, after updating dependencies:

```sh
npm --prefix web ci --ignore-scripts
npm --prefix web run sync:assets
```

Commit the updated bundles and licence notices. The browser build and this
command share the asset inventory in `web/scripts/assets.mjs`. CI runs
`npm --prefix web run check:assets` to detect drift without changing files.
Building the browser application does not rewrite these committed assets.
