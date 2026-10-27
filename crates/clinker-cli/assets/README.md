# Embedded viewer assets

`viewer/` is generated from `web/src/` and pinned D3/clustermap dependencies.
Cargo embeds this snapshot without running frontend tools. Edit the sources,
then refresh from the repository root:

```sh
npm --prefix web ci --ignore-scripts
npm --prefix web run sync:assets
npm --prefix web run check:assets
```

Commit the complete snapshot and licence notices, including renamed or removed
vendor chunks. The separate portable HTML exporter uses the bundled renderer
and D3. See [web development](../../../web/README.md).
