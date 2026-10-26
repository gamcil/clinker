#!/usr/bin/env sh
set -eu

# Preserve paths relative to the caller's working directory.
exec node "$(dirname "$0")/scripts/package.mjs" "$@"
