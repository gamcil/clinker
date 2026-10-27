# Python maintenance and PyPI releases

The Python implementation is maintained on `python-maintenance`. Its baseline
is the existing `v0.0.32` release, commit
`a126ee5306982858809981b79d6af3663cd2c56c`. The Python source and packaging on
`rust-port` matched that release when the maintenance branch was created.

Keep future Python fixes, dependency updates, package version changes and PyPI
releases on that branch. The Rust implementation is prepared on `rust-port` for
merging into `master`. The Rust branch no longer contains the Python
implementation, packaging files, Dockerfile, or PyPI workflow. Check out
`python-maintenance` before building or publishing a Python package.

## Publishing contract

The Python publishing workflow is triggered only by pushed `python-v*` tags.
Branch pushes and creating GitHub releases do not trigger this new workflow.
Keep the workflow on the maintenance branch when it is removed from the future
Rust-only branch.

Before building, the workflow requires:

- The tagged commit to be reachable from the remote `python-maintenance` branch.
- The tag to exactly equal `python-v` followed by `clinker.__version__`.
- A single literal, non-empty `__version__` string in `clinker/__init__.py`.

The source version is read without importing the package or installing its
runtime dependencies. The workflow then builds the wheel and source archive,
checks their metadata with Twine, and uploads them using the existing
`PYPI_TOKEN` repository secret. The token is supplied only to the upload step.

Python release tags use `python-vX.Y.Z`. Existing `v0.0.*` tags remain historical;
do not move or recreate them. Rust releases use `vX.Y.Z` from `master`, starting
with the planned `v0.1.0`. These tags do not match the maintenance branch's Python workflow.
Avoid creating new releases from historical commits containing the old broad
release-triggered workflow: workflow files are versioned with those commits.

## Future maintenance release

1. Check out `python-maintenance` and update it from the remote.
2. Make the Python fix and bump `__version__` to a new, unpublished version.
3. Review the change and check the Python package for the affected behaviour.
4. Build both distributions and check their metadata in an isolated environment:

   ```sh
   python -m pip install build twine
   python -m build
   python -m twine check --strict dist/*
   ```

5. Commit the fix and version bump, then push the maintenance branch before the
   tag. For example, for a future `0.0.33` release:

   ```sh
   git push origin python-maintenance
   git tag -a python-v0.0.33 -m "Python maintenance release 0.0.33"
   git push origin refs/tags/python-v0.0.33
   ```

Pushing the tag starts the PyPI upload; a GitHub release is optional and does
not act as a publishing approval. Confirm the Actions run and PyPI artifacts
before announcing the release. Failed validation must be corrected before
publishing. Existing PyPI versions cannot be replaced by rerunning the workflow.

Branch creation and this workflow migration do not bump the version, create a
release tag or upload a package. The existing PyPI secret and GitHub execution
still need to be confirmed during the first actual maintenance release.
