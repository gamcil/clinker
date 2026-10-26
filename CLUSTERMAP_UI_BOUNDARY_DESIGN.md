# Clustermap agent prompt: scoped UI styling and integration contract

## Task

Refactor clustermap's editor and other interactive controls to have a clear,
dependency-free styling boundary. Keep the editor as an optional, reusable
component distributed by clustermap. Use scoped CSS, explicit classes and a
documented theme contract. Do not add DaisyUI, Tailwind or another style library.

This is an implementation prompt. Inspect the current repository and its local
instructions before editing. The consumer inspected when writing this document
was Clinker, using clustermap commit
`a524445f5f656b684cf7fcd22c2286ad9d691bcd`; account for subsequent changes.

Companion prompt: `CLINKER_UI_BOUNDARY_DESIGN.md` in the Clinker repository.

## Ownership

| Layer | Responsibility |
| --- | --- |
| Clustermap core | Rendering, configuration, selection, ordering and annotation panel models |
| Clustermap UI/editor | Editor controls, selection toolbar, context menus, interaction, navigation and validation of figure settings |
| Host application | Analysis, file selection and parsing, cluster matching, project workflows, application navigation and UI appearance preference |
| Integration contract | Public methods/events, theme properties and host-controlled placement/sizing |

Clinker constructs annotation data from CSV/TSV; clustermap displays and edits
the resulting panels. Keep file import workflows in the host.

## Findings to address

- The editor already has shared styles in `src/editor/shared/controls.mjs` and
  theme handling in `src/editor/shared/theme.js`. Extend or simplify these
  rather than creating a competing control system.
- Editor components use nested shadow roots. This makes injecting host CSS an
  unsuitable integration mechanism, but does not prevent sharing internal CSS.
- The theme helper declares public `--cm-*` properties on each component host,
  with `--clinker-*` hooks and configuration-derived fallbacks. A declaration on
  a child host can mask an inherited public property; fix this inheritance trap.
- Clinker currently restyles the light DOM selection toolbar and reaches through
  editor shadow roots to navigate to tree settings. Both should become unnecessary.

## Styling architecture

1. Keep the editor's outer Shadow DOM boundary. Audit nested roots, retaining
   them where independent reuse or isolation warrants them. Do not flatten the
   entire editor merely to make host selectors reach its controls.
2. Consolidate shared controls and interaction states into a small internal CSS
   foundation. Component styles own layout and specialised presentation.
3. Use explicit classes for button, menu-item, field, section and state variants.
   Avoid broad document selectors, app-specific selectors and specificity fixes.
4. Scope light DOM UI under an owned root such as `.cm-ui`. Use that scope for
   floating menus, selection controls and other library-created HTML. Keep those
   rules away from the scientific SVG and unrelated host elements.
5. Ensure body-mounted/popover UI receives the same appearance as its owning
   chart. Document how inherited tokens reach it; do not assume a body-mounted
   element inherits variables from the chart container. Support two charts with
   different appearances on the same page.
6. Keep menu rows borderless, with hover, selected and visible keyboard-focus
   treatments. Input boundaries and container boundaries should be intentional.
7. Share stylesheet definitions across components; do not maintain parallel
   copies of the same button and field rules.

## Public appearance contract

Adopt a compact documented set of public properties. Proposed baseline:

| Properties | Purpose |
| --- | --- |
| `--cm-surface`, `--cm-subtle`, `--cm-hover` | UI surfaces and interaction backgrounds |
| `--cm-foreground`, `--cm-muted` | Normal and secondary text |
| `--cm-border`, `--cm-accent`, `--cm-focus-ring` | Boundaries, active controls and focus |
| `--cm-font` | UI typography |
| `--cm-control-height`, `--cm-control-radius` | Compact control dimensions |

Audit existing names before finalising this list. Reuse meaningful existing
names and add semantic tokens for error/success/disabled treatments where needed.
Publish the exact supported list and defaults in the repository documentation.

- Public overrides must inherit through nested components. Prefer consuming a
  public property with a private fallback, e.g.
  `color: var(--cm-foreground, var(--_cm-default-foreground));`, rather than
  assigning a fallback to `--cm-foreground` on every nested host.
- Standalone controls need complete, accessible defaults. Define a predictable
  precedence: explicit public override, explicit UI appearance mode, automatic
  system preference when no UI mode is supplied.
- UI appearance must be independently controllable from figure colours,
  annotation palettes, SVG background and exported figure configuration. Audit
  the current coupling to `config.plot.theme`; preserve existing figure settings
  and document compatibility behaviour rather than silently redefining them.
- Make colour-scheme explicit for native form controls and propagate the chosen
  UI mode consistently to nested and floating controls.
- Replace the app-specific `--clinker-*` contract with `--cm-*`. If compatibility
  aliases are necessary, keep them in one documented adapter with clear
  precedence; consumers must be able to migrate completely away from them.

## Public editor navigation

Add a supported way to open an appearance section, including tree settings.
An API such as `editor.openAppearance({ section: "tree", focus: true })` is a
proposal, not an existing API. Choose names consistent with the current editor.

The method must open the editor, select the correct tab, clear or adjust filters
that hide the target, reveal the target section, and focus its relevant control
when requested. Handle rendering timing within the editor. Existing annotation
and selection entry points should continue to work.

Document method arguments, readiness, unknown-section behaviour, return value
and open/close events. Keep events sufficient for the host to synchronise its
Edit figure button without inspecting editor markup. Reuse existing chart change
events rather than adding duplicate state ownership.

## Distribution and compatibility

- Preserve an optional editor entry point and independent core use. The core
  must not start loading editor CSS or registering editor elements implicitly.
- Ship editor styles with the editor; embedding applications need no CSS build
  tool or external stylesheet service.
- Preserve existing entry points and custom element names unless a separately
  justified migration is required. Renaming `<clinker-editor>` is outside scope.
- Update the build/package asset handling if stylesheet/module layout changes.
- No CDN requests, new runtime styling dependencies or separate design-system
  package are required.

## Acceptance and handoff

Document evidence for the following in standalone and embedded examples:

- Light/dark UI, native inputs, disabled states, focus and selected states are legible.
- Tokens set on a host wrapper reach nested editor controls and floating UI.
- Two chart instances can use different UI palettes without leakage.
- Ordinary host button/input rules do not alter shadow editor internals; library
  light DOM rules do not alter unrelated host controls.
- Batch selection toolbar and menus need no Clinker component-internal overrides.
- Tree navigation works from a closed editor and with a hiding filter active.
- UI appearance changes leave figure configuration, saved data and exports intact.
- Renderer-only and packaged editor consumption remain functional.

Commit at clean boundaries: internal scoped styles; public appearance contract
and compatibility; navigation API and documentation. Report the resulting commit
hash, exact API/token contract, migration notes and asset changes to the Clinker
agent. Do not edit a consumer's installed `node_modules` as the source of truth.
