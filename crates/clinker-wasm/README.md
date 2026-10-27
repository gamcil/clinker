# clinker-wasm

Internal WebAssembly bindings for the clinker browser application. This workspace
crate is not published to crates.io. Build the browser application from the
repository root with `sh web/build.sh`; the resulting static app is in `web/dist/`.

The analysis implementation lives in `clinker-core`. Native CLI users do not
need this crate or a WASM toolchain.
