# Rust release procedure

The hosted application is the primary interface. The same workflow validates
pushes and pull requests, deploys `rust-port` to Pages, and builds release
archives for tags matching `rust-v*`.

## Before tagging

- Set the workspace version in `Cargo.toml` and the `clinker-core` dependency
  version in `crates/clinker-cli/Cargo.toml`; refresh `Cargo.lock` if needed.
- Review [FEATURES.md](FEATURES.md), [RELEASE_NOTES.md](RELEASE_NOTES.md), and installation
  instructions in [CLI.md](CLI.md) and [web/DEPLOYMENT.md](web/DEPLOYMENT.md).
- Confirm the checks and hosted Pages build succeed for the intended commit.
- Manually check the hosted app in current Safari. Automated upload, analysis,
  cancellation, editing, project reopening and export checks run in Chromium,
  Firefox and WebKit; WebKit coverage does not replace testing actual Safari.
- Confirm native tests and release-binary smoke checks pass on Windows, Linux
  and Apple Silicon macOS. Tag builds package those binaries; downloaded
  archives still need extraction checks on their target operating systems.

The deployment branch remains `rust-port`. Moving the final application to
`master` requires updating the two Pages conditions in the workflow and the
branch name in this document. Do not change that routing just to make a tag.

## Publish

Create and push `rust-vX.Y.Z` on the intended commit, with `X.Y.Z` exactly
matching the workspace version. A tag push publishes a GitHub release after
native tests, lint, WASM tests, and packaged browser tests pass. It does not
deploy Pages. The release contains three native archives, `clinker-web.tar.gz`,
and `SHA256SUMS`. Rust releases are excluded from Python package publishing.

Python maintenance is isolated on `python-maintenance`, starting from the
Python `v0.0.32` release. Its publishing workflow accepts only pushed
`python-v*` tags from that branch's history, with a matching Python package
version. It does not subscribe to GitHub release events. See
[Python maintenance releases](../docs/python-maintenance.md). Keep the existing
`rust-v*` tag convention until the later Rust-first release cutover.

Download the artifacts, verify their checksums, and check that the static web
archive works when served over HTTP. GitHub execution, cross-platform archive
extraction and the actual public Pages deployment must be confirmed on GitHub;
local checks do not substitute for those steps.

