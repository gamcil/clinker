# Browser app

Static GenBank and GFF3/FASTA comparison using the shared Rust core through
WASM and Web Workers. Input files stay in the browser. Clustermap supplies the
renderer and editor; no backend service is required.

## Build and preview

Requires Node/npm, rustup and wasm-pack. From the repository root:

```sh
rustup target add wasm32-unknown-unknown
cargo install wasm-pack --locked
sh web/build.sh
python3 -m http.server --directory web/dist
```

Open the reported HTTP URL. `build.sh` installs pinned npm dependencies and
uses rustup's toolchain. Build output is replaced on each run.

```sh
sh web/package.sh /tmp/clinker-web
```

This exports the built application for hosting. See [deployment](DEPLOYMENT.md)
and [features and data imports](../docs/features.md).

## Tests

```sh
npm --prefix web test
npm --prefix web run test:wasm
cd web
npx playwright install chromium
npm run test:browser
```

WASM and browser tests require a completed build. Set `CLINKER_TEST_BROWSER` to
`firefox` or `webkit` to use another installed engine. CI runs Chromium on normal
pushes and all three engines for manual runs and release tags.

## Sources and assets

- `src/app.js`: browser entry and analysis orchestration.
- `src/local.js`: native viewer entry; uses the CLI's project/similarity endpoints.
- `src/viewer/`: shared figure lifecycle, controls, project validation and tree helpers.
- `src/tree-*.mjs`, `src/annotation-import.mjs`: import and matching workflows.
- `scripts/`: build assembly, packaging and the shared vendor inventory.
- `.build/`: intermediate WASM and native viewer output; `dist/`: static app.

Builds preserve upstream vendor paths and add a content-derived version to app
asset URLs. `build.json` records the packaged file inventory and hashes.
After changing viewer source or dependencies, run `npm --prefix web run
sync:assets` and commit the generated CLI snapshot. `check:assets` verifies it.
Cargo builds use that snapshot and need no Node.js. See
[architecture](../docs/architecture.md) for the UI integration boundary.
