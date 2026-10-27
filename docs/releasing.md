# Rust release procedure

The hosted application is the primary interface. The same workflow validates
pushes and pull requests, deploys `rust-port` to Pages, and builds release
archives for tags matching `v[0-9]*`.

## Before tagging

- Set the workspace version in the root `Cargo.toml` and the `clinker-core` dependency
  version in `crates/clinker-cli/Cargo.toml`; refresh `Cargo.lock` if needed.
- Run the Cargo source-package checks from a clean checkout. They test and
  install the extracted archives, including project/HTML exports and the embedded
  viewer. See [packaging](packaging.md#source-package-verification).
- Review [features.md](features.md), [release-notes.md](release-notes.md), and installation
  instructions in [cli.md](cli.md) and [web/DEPLOYMENT.md](../web/DEPLOYMENT.md).
- Confirm the checks and hosted Pages build succeed for the intended commit.
- Manually check the hosted app in current Safari. Automated upload, analysis,
  cancellation, editing, project reopening and export checks run in Chromium on
  ordinary pushes. Run the workflow manually to check Firefox and WebKit before
  tagging; release tags also run all three. WebKit does not replace actual Safari.
- Confirm native tests and archive extraction checks pass on Windows, Linux
  and Apple Silicon macOS. Every native job builds an optimized binary, packages
  it, extracts the archive and runs that extracted executable. Check a downloaded
  release on each supported OS as well.

## Repository cutover

Rust development and Pages deployment remain on `rust-port` while the code base
is finalized. The merge into `master` is deferred. The executable rename and
packaging work can be used independently of that merge.

When ready to make Rust the main repository version:

1. Push the prepared Rust branch and confirm all native, source-package and web
   checks pass on GitHub.
2. Change both Pages conditions in `.github/workflows/deploy-pages.yml` from
   `refs/heads/rust-port` to `refs/heads/master`, then merge into `master`. Keep
   `python-maintenance` as the Python/PyPI branch.
3. Confirm GitHub's default branch is `master`, Pages uses GitHub Actions, and
   the `github-pages` environment permits deployment from `master`.
4. Confirm the merged commit passes CI and deploys the intended hosted app.
5. Review the first Rust release at the existing workspace version `0.1.0`;
   create `v0.1.0` only when ready to publish. Set a different version before
   tagging if required.

Preparing this commit does not push a branch, merge, alter repository settings,
create a tag or publish a release. Those actions remain the actual cutover.

## Publish

Create and push `vX.Y.Z` on a commit already merged into `master`, with `X.Y.Z`
exactly matching the workspace version. CI checks both version and ancestry.
A tag push publishes a GitHub release after native tests, lint, source-package
installation checks, WASM tests, and packaged
browser tests pass. It does not deploy Pages. The release contains three native
archives, `clinker-web.tar.gz`, the core/CLI `.crate` archives, and `SHA256SUMS`.
Rust releases are excluded from Python package publishing.

The workflow does not publish to crates.io. Follow the separate core-then-CLI
[registry publication procedure](packaging.md#registry-publication) when ready.
The CLI package installs one executable named `clinker`. Native archives and
Cargo installations use the same executable name.

Python maintenance is isolated on `python-maintenance`, starting from the
Python `v0.0.32` release. Its publishing workflow accepts only pushed
`python-v*` tags from that branch's history, with a matching Python package
version. It does not subscribe to GitHub release events. See
[Python maintenance releases](python-maintenance.md). Existing Python `v0.0.*`
tags remain historical; do not move or recreate them.

Download the artifacts, verify their checksums, and check that the static web
archive works when served over HTTP. GitHub execution, cross-platform archive
extraction and the actual public Pages deployment must be confirmed on GitHub;
local checks do not substitute for those steps.

