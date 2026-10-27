# Saved projects and analysis sessions

The browser and native CLI use the `clinker-project` version **2** envelope for
editable figures with their analysis data. This retains the existing browser
format; no migration or version bump is needed for saved version 2 projects.

```json
{
  "format": "clinker-project",
  "version": 2,
  "analysis": {
    "format": "clinker-analysis",
    "version": 1,
    "identityCutoff": 0.3,
    "analysis": { "clusters": [], "links": [] },
    "groupDefinitions": []
  },
  "data": { "clusters": [], "links": [], "groups": [] },
  "config": {}
}
```

The empty arrays above illustrate the envelope; actual analysis records contain
parsed loci, genes and translations. A project contains:

| Field | Owner | Meaning |
| --- | --- | --- |
| `format`, `version` | Clinker | Project format discriminator |
| `analysis` | clinker-core | Validated version 1 analysis session, including translations |
| `data` | clustermap | Current editable figure, including annotations and any active tree |
| `config` | clustermap | Figure appearance settings |
| `state` (optional) | clustermap | Version 1 layout state, including cluster order and camera |
| `inactiveTree` (optional) | Clinker | Imported tree retained while another row order is active |
| `clusteringDendrogram` (optional) | Clinker | Retained clustering tree |
| `treeSource` (optional) | Clinker | `clustering` when the active tree is the computed dendrogram |

Interface System/Light/Dark preference belongs to the host application and is
not a figure setting. Loading a project does not change that preference.

## Preserving edits

The figure is the authority for visible names, geometry, annotations, tree edits,
ordering and deletions. Analysis retains the original parsed records and links.
Cluster, locus and gene IDs encode their original analysis indices; do not
renumber them after edits or reconstruct a saved figure from analysis on load.
Links and groups may retain references to genes deleted from the visible figure.

The shared JavaScript project module validates figure IDs against analysis,
coordinates, trees and layout, and retains unknown project extension fields on
export. The Rust `ClinkerProject` model validates the envelope and analysis
session and preserves figure data, config, layout and extra fields as JSON.
Renderer-specific validation remains in the shared JavaScript viewer. Native
code must not deserialize those fields into a narrower schema that drops edits.

## CLI and browser

```sh
clinker-rs input-a.gbk input-b.gbk --project comparison.json
clinker-rs --session analysis.json --project comparison.json
```

Open `comparison.json` through **File → Open project** in the browser app.
It opens the saved figure and uses the analysis payload to recover similarity
information; it does not rerun protein alignments or replace the figure layout.
New CLI projects use default figure appearance and omit saved camera state.

`--session` continues to read/write `clinker-analysis` version 1 for native
analysis continuation. `--project` writes an editable project; it does not change
the session persistence format. Standalone `--plot *.json` remains raw renderer
data and is not an editable project. Portable `--plot *.html` uses the separate renderer-only native HTML exporter.
Bare `--plot` opens the shared local viewer, and `--open-project comparison.json`
restores a saved project there without redoing analysis.

`ClinkerProject::from_json` and `to_json` provide the native project reader and
writer for the native local viewer without recomputing an edited figure.
