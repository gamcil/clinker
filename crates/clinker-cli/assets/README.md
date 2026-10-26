# Embedded plot assets

The CLI embeds `d3.min.js` and `clustermap.min.js` to make `--plot *.html`
self-contained. These generated bundles are deliberately vendored with the
crate so native Cargo builds do not depend on Node, npm, or network access.

Their source versions are pinned by `web/package-lock.json`. After
updating either dependency, run `npm ci --ignore-scripts` in `web`, copy
the corresponding distribution bundles into this directory, and update the
license copies when necessary. Rust CI compares these files with the web build
to prevent the two distributions from drifting.
