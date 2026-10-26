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

Resolved upstream in refactor commit
`f75cf112b7aee0c6315e1054c12645705578f6d9`. Clinker now pins that commit and
uses its published browser/editor assets and native CLI bundle.

The patch above is retained for historical context; do not apply it again.
