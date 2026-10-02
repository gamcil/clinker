#!/usr/bin/env sh
set -eu

cd "$(dirname "$0")"

# package-lock.json pins the refactored clustermap.js commit and D3 version.
npm ci --ignore-scripts

rustup_cargo=$(rustup which cargo)
rustup_bin=$(dirname "$rustup_cargo")
wasm_pack="${CARGO_HOME:-$HOME/.cargo}/bin/wasm-pack"

if [ ! -x "$wasm_pack" ]; then
  wasm_pack=$(command -v wasm-pack || true)
fi
if [ -z "${wasm_pack:-}" ]; then
  echo "wasm-pack was not found; install it with: cargo install wasm-pack" >&2
  exit 1
fi

# Prefer rustup's cargo, which has the wasm32 target, over Homebrew's cargo.
wasm_pack_bin=$(dirname "$wasm_pack")
PATH="$wasm_pack_bin:$rustup_bin:$PATH"
"$wasm_pack" build ../crates/clinker-wasm --target web --out-dir ../../web/pkg --out-name clinker_wasm
cp node_modules/d3/dist/d3.min.js .
cp node_modules/clustermap.js/dist/clustermap.min.js .
cp node_modules/clustermap.js/dist/editor.mjs clustermap-editor.mjs
cp node_modules/clustermap.js/dist/editor.css clustermap-editor.css
mkdir -p examples
cp ../../examples/*.gbk examples/
