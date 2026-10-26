# Clustermap agent: annotation legends overlap neighbouring panels

## Resolution

Implemented upstream in `8771037325631866d46316bd3379233677283c31` and imported
into Clinker through refactor commit `4a957c772f2e8e062f8c85ab21fb22b22811c095`.
The layout now reserves complete measured footprints on both sides, uses those
bounds for hover/figure geometry, and invalidates cached font measurements.
The original reproduction and requested acceptance scenarios below are retained
as historical context. Browser verification in Clinker has not yet been run.

## Reproduction

Clinker is pinned to `a524445f5f656b684cf7fcd22c2286ad9d691bcd`.
Plot the five example clusters and import `examples/annotations.tsv` as separate
tracks. The categorical organism legend extends underneath the next categorical
legend and the numeric bar-chart axis. The same issue can affect long headings
and panel titles. It also occurs in the current uncommitted local clustermap
source inspected on 2026-10-26.

## Cause

In `src/annotations/model.mjs`, `layoutAnnotationPanels()` computes
`panel.extentWidth` from header/footer text bounds, but advances the next panel
with `x += panel.width`. The reuse path also advances by `panel.width`.
Thus figure/export bounds include text that inter-panel spacing ignores.

## Requested fix in the library

Make automatic placement reserve the complete rendered horizontal footprint of
each panel, including cells/bars, keys, title, rotated headings and legend text.
Use measured text bounds, including left overhang where applicable. Preserve the
meaning of cell/bar width: increasing spacing must not stretch data marks.

Advancing by `extentWidth` in both the normal and cached paths is the minimal
correction for right overhang. Account for left overhang as well in the full
solution. Ensure cached geometry is invalidated when font measurements, fonts,
labels or presentation change. Keep explicit gap as additional separation.
Do not use a fixed large gap or hide legends as the default workaround.

Carry the corrected positions through figure bounds, SVG/Canvas/WebGPU rendering,
exports, minimap, hover hit areas, drag previews and panel reordering. If some
existing intentional overlapping/manual placement mode exists, document how it
differs from automatic layout.

## Acceptance scenarios

- Separate organism/class/length/count/yes-no tracks from the example table.
- Long categorical labels, long titles, rotated column headings, and bar axes.
- Changing a category label, font size or font; hide/show; reorder; undo/redo.
- Repeat the same layout using cached panels, then change a preceding panel.
- Pan/zoom and verify hover frames correspond to the new panel positions.
- Export and renderer switching preserve the same non-overlapping layout.

## Separate Clinker integration fix

Missing annotation hover/edit UI had a separate cause: Clinker's build and static
packaging omitted `dist/annotationOverlay.mjs`, which the renderer dynamically
imports. Clinker now copies and packages this module alongside the renderer and
existing `editor/` chunks, and handles `clustermap:edit-annotation` using the public
`editor.editAnnotations(panelUid)` entry point.

The library layout source has not been edited from Clinker. The local clustermap
working tree contains extensive independent changes. Please implement the layout
correction there and return a committed version for Clinker to pin and integrate.
