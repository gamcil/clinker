# clinker-wasm

Internal WebAssembly bindings for the clinker browser application. This workspace
crate is not published to crates.io. Build the browser application from the
repository root with `sh web/build.sh`; the resulting static app is in `web/dist/`.

The analysis implementation lives in `clinker-core`. Native CLI users do not
need this crate or a WASM toolchain.

Bindings reuse the core `GeneRef`, `Link` and `GeneGroupDefinition` types.
Browser-specific structures describe owned file bytes, parsing controls, metadata
without translations, protein worker inputs and JavaScript message envelopes.
Separating metadata from proteins keeps sequences out of the layout payload
sent between workers.
