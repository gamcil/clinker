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

Before the first registry publication, the CLI needs a local override for its
unpublished core dependency. These are the same steps used by CI, from a clean
checkout (Bash):

```sh
package_dir=$(mktemp -d)
cargo package --locked --no-verify -p clinker-core -p clinker-cli --target-dir "$package_dir"
for archive in "$package_dir"/package/*.crate; do
  tar -xzf "$archive" -C "$package_dir"
done
core=("$package_dir"/clinker-core-*)
cli=("$package_dir"/clinker-cli-*)
patch="patch.crates-io.clinker-core.path=\"${core[0]}\""
cargo test --manifest-path "${cli[0]}/Cargo.toml" --config "$patch" --target-dir "$package_dir/target"
cargo install --locked --debug --path "${cli[0]}" --config "$patch" --target-dir "$package_dir/target" --root "$package_dir/install"
(cd "$package_dir" && install/bin/clinker --version)
```

Cargo updates the extracted CLI lockfile for the core path override during the
test step; installation then uses that lockfile with `--locked`. Builds use a
fresh temporary directory and leave the checkout unchanged.

`--no-verify` skips Cargo's automatic registry verification of the unpublished
peer; the extracted-package tests perform verification instead. After the core
version is available on crates.io, replace these temporary extraction/override
steps with `cargo package --locked -p clinker-core -p clinker-cli`.

CI uploads the `.crate` files from its temporary package directory for tagged
GitHub releases. Debug installation keeps the source check fast; native jobs
separately build optimized binaries on each supported OS.

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
