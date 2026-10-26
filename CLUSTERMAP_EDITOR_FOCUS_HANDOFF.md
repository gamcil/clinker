# Clustermap editor focus handoff

## Problem

In Firefox, opening tree appearance from Clinker's Row order menu focuses the
editor, then closing the menu exposes the raster canvas beneath the pointer.
The canvas pointer-enter handler calls `focus()` and takes focus back from the
editor. This breaks keyboard navigation and the tree settings browser check.

The window Space handler also sees a shadow host as `event.target`, so it can
incorrectly start plot panning while a control inside the editor handles Space.

## Proposed upstream patch

`CLUSTERMAP_EDITOR_FOCUS.patch` contains only the raster interaction source change
and its unit regression. It was prepared against refactor commit `4218e0a`.

From the clustermap repository:

```sh
git apply --check /Users/gamcil/repos/clinker/CLUSTERMAP_EDITOR_FOCUS.patch
git apply /Users/gamcil/repos/clinker/CLUSTERMAP_EDITOR_FOCUS.patch
```

Review/adapt the patch against the current source before applying if it has moved.
Regenerate distribution assets using clustermap's normal build process.

Expected behavior:

- Hover may focus an idle canvas, but preserves focus in editor and other controls.
- A primary pointer press explicitly focuses the canvas.
- Space inside a control, including a shadow DOM control, remains with that control.

Verify the raster unit regression and the Firefox flow: open Row order → Tree
appearance, then use the focused tree control by keyboard. Check that canvas
clicking and Space panning still work when the canvas has focus.

## Integration status

An earlier attempt already created upstream branch `fix/preserve-editor-focus`
at `c30feeaba3dfaac03e10a328218dd4b183db08dc`. It includes this source patch and
rebuilt distribution files. It has not been merged into `refactor`; use either
that commit or this patch, without applying both.

Clinker's integration commit `21f5164` is reverted in the accompanying change,
restoring its previous pin `4a957c772f2e8e062f8c85ab21fb22b22811c095` and vendored
bundle. After the fix is integrated upstream, update Clinker's dependency and
generated assets from the resulting refactor commit. Until then the Firefox
focus regression remains unresolved in that pinned version.
