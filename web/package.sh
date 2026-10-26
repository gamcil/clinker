#!/usr/bin/env sh
set -eu

destination=${1:?Specify a destination for the static application}
case "$destination" in
  /*) ;;
  *) destination="$(pwd)/$destination" ;;
esac
cd "$(dirname "$0")"
mkdir -p "$destination"
cp index.html app.js worker.js style.css input-files.mjs alignment.mjs workflow.mjs project.mjs tree.mjs tree-import.mjs tree-input.mjs annotation-import.mjs \
  d3.min.js clustermap.min.js clustermap-editor.mjs annotationOverlay.mjs "$destination/"
cp -R editor pkg examples "$destination/"
cp ../LICENSE ../crates/clinker-cli/assets/LICENSE.* "$destination/"
cp DEPLOYMENT.md "$destination/README.md"
