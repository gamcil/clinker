# Regression testing

Tests protect supported workflows and known failures. There is no coverage
percentage or test-count target. Prefer representative fixtures and observable
results; keep focused algorithm tests where a workflow would hide the failure.

## Ownership

| Owner | What to verify |
| --- | --- |
| Rust core | Parsing and translation, alignment, filtering, grouping, ordering, layout and saved-analysis invariants |
| CLI | Input discovery, session continuation, exports and preservation of existing files after failure |
| WASM adapter | Conversion of browser inputs/options and results across the JavaScript boundary |
| Clustermap | Renderer conformance, chart interactions, editor behavior, annotations, tree edits and undo/redo |
| Clinker browser app | Analysis/cancellation, import and matching, editor integration, project save/reopen and exports |
| Distribution | Installation from source archives and working native/web release artifacts |

An algorithm's detailed cases belong to its owner. Adapters need representative
boundary checks, not a second implementation of the algorithm's test suite.
Clinker exercises an edit and project round-trip against its pinned clustermap;
exhaustive renderer, dragging and tree-operation checks belong upstream.

## Keeping the suite small

- Add a regression for a reproducible bug or a supported behavior that lacks a
  useful check. Extend an existing scenario when the setup and contract match.
- Share ordinary setup, not a general-purpose test framework. Keep assertions
  close to the actions and avoid source-code patterns, private object shapes,
  method inventories and precise cosmetic properties as test contracts.
- Keep independent workflows separately diagnosable. Combining unrelated cases
  into one large test solely to reduce the count is not simplification.
- Run Chromium browser workflows on ordinary pushes and pull requests. Release
  tags and manual workflow runs also exercise Firefox and WebKit. Native builds
  and filesystem regressions continue on Linux, macOS and Windows.
