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

## Source package verification

CI uses Cargo to create the source archives, then tests and installs the CLI
from their extracted contents outside the checkout. The CLI workflow regression
covers project/HTML exports and serving the embedded viewer and saved project.
Native jobs also extract their platform archive and run its executable outside
the checkout.

Now that `clinker-core` and `clinker-cli` are published, Cargo can resolve the
CLI's core dependency from the registry. From a clean checkout, package and
verify both published crates with:

```sh
cargo package --locked -p clinker-core -p clinker-cli
```

Cargo creates the `.crate` archives under `target/package/` and verifies each
archive builds from its extracted contents. The CI source-package workflow also
tests and installs the extracted CLI, including its embedded viewer.

CI uploads the `.crate` files from its temporary package directory for tagged
GitHub releases. Debug installation keeps the source check fast; native jobs
separately build optimized binaries on each supported OS.

## Registry publication

Registry publication is a separate release action, not part of the GitHub
workflow. Version `0.1.0` of both public crates has been published. For each
subsequent release, update the workspace version and the CLI's `clinker-core`
dependency version together, then publish from a clean release commit in order:

```sh
cargo publish --dry-run --locked -p clinker-core
cargo publish --locked -p clinker-core
cargo info clinker-core
cargo publish --dry-run --locked -p clinker-cli
cargo publish --locked -p clinker-cli
```

Wait for the core version to appear in the registry before publishing the CLI.
Cargo polls the registry after uploading; see
[Cargo publish](https://doc.rust-lang.org/cargo/commands/cargo-publish.html).
Then confirm a registry installation with `cargo install --locked clinker-cli`.
