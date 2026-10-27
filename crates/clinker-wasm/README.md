# clinker-wasm

Internal WebAssembly adapter for `clinker-core`, used by the browser app and
not published to crates.io. Build from the repository root with `sh web/build.sh`.

The adapter converts browser inputs and worker messages; analysis stays in the
shared core. Native CLI builds do not require WASM tooling.
