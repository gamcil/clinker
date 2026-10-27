# Rust distribution packages

## Outputs

| Output | Contents | Consumer requirements |
| --- | --- | --- |
| `clinker-core-VERSION.crate` | Library source, README, licence, examples and package-local regression fixtures | Rust and a native build toolchain |
| `clinker-cli-VERSION.crate` | CLI source, Rust build script, complete embedded viewer, tests, README and licences | Rust and a native build toolchain |
| Native `.tar.gz` / `.zip` | One executable, CLI/project-format documentation and licence notices | Supported OS and a current browser for the viewer |
| `clinker-web.tar.gz` | Tested static browser app, WASM, vendor modules, examples and licences | An HTTP static host and a supported browser |

`clinker-wasm` has `publish = false`; its bindings are distributed through the
static application. The CLI package is named `clinker-cli` and installs one executable, `clinker`.
The workflow's `CLINKER_BINARY` value is checked against Cargo metadata to keep
archive packaging consistent.

## Cargo contents

Each published crate has explicit include rules and its own README/licence.
Tests use crate-local fixtures, including frozen parser examples; they do not
reach into a sibling crate or the repository-level demonstration inputs.

The CLI's committed `assets/viewer/` contains its complete viewer. `build.rs`
indexes those files and emits `include_bytes!` references into Cargo's output
directory. It never modifies source files or invokes frontend tools. Both Cargo
and native binary consumers can build/use the viewer without a web checkout.

The core dependency declares both `path` and `version`. Cargo removes the path
when packaging for a registry, so publish the core before publishing the CLI.
See [Cargo packaging](https://doc.rust-lang.org/cargo/commands/cargo-package.html).

## Isolated package verification

From a clean checkout:

```sh
node scripts/check-packages.mjs
```

For local development, `--offline` uses cached Cargo dependencies and
`--allow-dirty` explicitly permits packaging modified source. CI uses neither.
The script requires Node, Cargo and `tar` as maintainer tools; those tools do not
become dependencies of an installed application.

The check:

1. Packages core and CLI together with the workspace lockfile and explicit file lists.
2. Extracts both archives into a temporary workspace outside the checkout.
3. Overrides only the unpublished core dependency to the extracted peer archive,
   keeping every third-party locked version. It never uses the checkout's core.
4. Builds all packaged targets and performs a debug `cargo install --locked` into
   a temporary installation root.
5. Uses the installed executable to generate a project and portable HTML.
6. Removes both extracted source trees, starts the installed viewer, and checks
   the served project and the hash of every embedded asset.

Package creation uses `--no-verify` because verification is explicitly performed
on the extracted pair. This avoids requiring the first core version to exist in
the public registry and works around Cargo 1.94's temporary-registry checksum
failure for interdependent package verification. Debug installation keeps this
check fast; native release jobs separately build optimized binaries on each OS.

Accepted archives remain in `target/package-check/package/`. CI gates releases
and Pages on this check and adds the `.crate` files to tagged GitHub releases.
Native jobs also extract their platform archive and run its executable outside
the checkout on every run.

## Registry publication

Registry publication is a separate release action, not part of the GitHub
workflow. Confirm crate-name ownership and registry permissions before the first
publication. From the clean release commit, publish in order:

```sh
cargo publish --locked -p clinker-core
cargo publish --locked -p clinker-cli
```

Wait for the core version to become available before publishing the CLI. Cargo
polls the registry after uploading; see
[Cargo publish](https://doc.rust-lang.org/cargo/commands/cargo-publish.html).
Then confirm a registry installation with `cargo install --locked clinker-cli`.
The isolated archive check exercises installation before that first publication;
it does not claim to verify registry ownership or upload credentials.
