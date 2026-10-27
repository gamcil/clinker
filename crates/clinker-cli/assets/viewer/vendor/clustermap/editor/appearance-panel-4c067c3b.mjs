import { i as color, m as mountComponent, c as createTemplate, b as uiAppearanceStyles, l as lightControlStyles, f as defineSearchField, g as createActionRunner, u as upgradeProperty, e as colourHex } from './search-field-595834d9.mjs';
import { a as applyEditorTheme, e as editorThemeStyles } from './theme-8bd47f2f.mjs';

const TREE_COLOURING_DEFAULTS = {
  mode: "none",
  cutHeight: "auto",
  palette: ["#267c9e", "#d7752d", "#5c8f39", "#9963a5", "#c34d61", "#9e8433"],
  showCutLine: false,
};
const metricsCache = new WeakMap();
const stylesCache = new WeakMap();
const fail = (message) => { throw new TypeError(`Invalid tree presentation: ${message}`); };
const object = (value) => value && typeof value === "object" && !Array.isArray(value);
const own = (value, key) => Object.prototype.hasOwnProperty.call(value, key);
function fields(value, allowed, where) {
  if (!object(value)) fail(`${where} must be an object.`);
  for (const key of Object.keys(value)) if (!allowed.includes(key)) fail(`${where}.${key} is unsupported.`);
}
function colour(value, where) {
  if (typeof value !== "string" || !color(value)) fail(`${where} must be a valid colour.`);
}
function positive(value, where) {
  if (!Number.isFinite(value) || value <= 0) fail(`${where} must be finite and positive.`);
}
function validateTreeColouring(policy) {
  if (policy === undefined) return;
  fields(policy, ["mode", "cutHeight", "palette", "showCutLine"], "colouring");
  if (own(policy, "mode") && !["none", "cutHeight"].includes(policy.mode)) fail("colouring.mode must be none or cutHeight.");
  if (own(policy, "cutHeight") && policy.cutHeight !== "auto" && (!Number.isFinite(policy.cutHeight) || policy.cutHeight < 0)) {
    fail('colouring.cutHeight must be "auto" or a finite non-negative number.');
  }
  if (own(policy, "showCutLine") && typeof policy.showCutLine !== "boolean") fail("colouring.showCutLine must be boolean.");
  if (own(policy, "palette")) {
    if (!Array.isArray(policy.palette) || !policy.palette.length) fail("colouring.palette must be a non-empty colour array.");
    policy.palette.forEach((entry) => colour(entry, "colouring.palette"));
  }
}

/** Biological distances and canonical ID order survive child rotation unchanged. */
function treeMetrics(tree) {
  if (metricsCache.has(tree)) return metricsCache.get(tree);
  const nodes = new Map(tree.nodes.map((node) => [node.uid, node]));
  const parents = new Map();
  const distances = new Map([[tree.root, 0]]), depths = new Map([[tree.root, 0]]);
  let completeLengths = true, nonNegative = true;
  const visit = (uid) => {
    for (const child of nodes.get(uid).children || []) {
      const length = nodes.get(child).branchLength;
      const finite = Number.isFinite(length);
      if (!finite) completeLengths = false;
      if (length < 0) nonNegative = false;
      parents.set(child, uid);
      distances.set(child, distances.get(uid) + (finite ? length : 0));
      depths.set(child, depths.get(uid) + 1);
      visit(child);
    }
  };
  visit(tree.root);
  const tips = tree.nodes.filter((node) => !node.children);
  const height = Math.max(0, ...tips.map((node) => distances.get(node.uid)));
  // Absolute 1e-12 plus relative 1e-9 tolerance, measured in input distance units.
  const tolerance = 1e-12 + 1e-9 * height;
  const ultrametric = completeLengths && nonNegative && tips.every((node) => Math.abs(distances.get(node.uid) - height) <= tolerance);
  const result = {
    nodes, parents, distances, depths, maxDepth: Math.max(...depths.values()), height, completeLengths, ultrametric,
    canonical: [...tree.nodes].sort((a, b) => a.uid < b.uid ? -1 : a.uid > b.uid ? 1 : 0),
  };
  metricsCache.set(tree, result);
  return result;
}

function validateBranch(style, where) {
  fields(style, ["colour", "width"], where);
  if (own(style, "colour")) colour(style.colour, `${where}.colour`);
  if (own(style, "width")) positive(style.width, `${where}.width`);
}

/** Validate optional data fields and the effective tree/config colouring policy. */
function validateTreePresentation(tree, config = {}) {
  validateTreeColouring(config.colouring);
  if (!tree) return;
  if (tree.distance !== undefined) {
    fields(tree.distance, ["kind", "label"], "distance");
    if (tree.distance.kind !== "mergeHeight") fail("distance.kind must be mergeHeight.");
    if (own(tree.distance, "label") && typeof tree.distance.label !== "string") fail("distance.label must be text.");
  }
  const presentation = tree.presentation;
  if (presentation !== undefined) {
    fields(presentation, ["colouring", "subtrees", "nodes"], "presentation");
    validateTreeColouring(presentation.colouring);
    const nodes = treeMetrics(tree).nodes;
    for (const scope of ["subtrees", "nodes"]) {
      if (presentation[scope] === undefined) continue;
      if (!object(presentation[scope])) fail(`presentation.${scope} must be a node-ID dictionary.`);
      for (const [uid, record] of Object.entries(presentation[scope])) {
        const node = nodes.get(uid);
        if (!node) fail(`presentation.${scope} refers to unknown node '${uid}'.`);
        fields(record, scope === "subtrees" ? ["branch"] : ["branch", "label"], `${scope}.${uid}`);
        if (scope === "subtrees" && !node.children) fail(`subtree '${uid}' must be an internal node.`);
        if (record.branch !== undefined) {
          if (scope === "nodes" && uid === tree.root) fail("the display root has no incoming edge to style.");
          validateBranch(record.branch, `${scope}.${uid}.branch`);
        }
        if (record.label !== undefined) {
          if (!node.children) fail(`annotation '${uid}' must belong to an internal node.`);
          const label = record.label;
          fields(label, ["text", "colour", "fontSize", "fontWeight", "fontStyle"], `nodes.${uid}.label`);
          if (own(label, "text") && typeof label.text !== "string") fail(`annotation '${uid}' text must be a string.`);
          if (own(label, "colour")) colour(label.colour, `annotation '${uid}' colour`);
          if (own(label, "fontSize")) positive(label.fontSize, `annotation '${uid}' fontSize`);
          if (own(label, "fontWeight") && ![400, 700].includes(label.fontWeight)) fail("annotation fontWeight must be 400 or 700.");
          if (own(label, "fontStyle") && !["normal", "italic"].includes(label.fontStyle)) fail("annotation fontStyle must be normal or italic.");
        }
      }
    }
  }
  const policy = { ...TREE_COLOURING_DEFAULTS, ...config.colouring, ...presentation?.colouring };
  if (policy.mode === "cutHeight") {
    if (tree.distance?.kind !== "mergeHeight") fail("cutHeight colouring requires distance.kind = mergeHeight.");
    if (!treeMetrics(tree).ultrametric) fail("cutHeight colouring requires complete lengths and equal root-to-tip distances (absolute 1e-12 + relative 1e-9 tolerance).");
  }
}

/** Resolve properties once per presentation/config revision, never per tween. */
function resolveTreePresentation(tree, config = {}) {
  const defaults = { ...TREE_COLOURING_DEFAULTS, ...config.colouring };
  const signature = JSON.stringify([config.branch?.colour, config.branch?.width, config.textColour, defaults]);
  const cached = stylesCache.get(tree);
  if (cached && cached.presentation === tree.presentation && cached.signature === signature) return cached.value;
  validateTreePresentation(tree, config);
  const metrics = treeMetrics(tree);
  const policy = { ...defaults, ...tree.presentation?.colouring };
  const threshold = policy.cutHeight === "auto" ? metrics.height * 0.7 : policy.cutHeight;
  const groups = new Map();
  if (policy.mode === "cutHeight" && threshold > 0) {
    const heightFor = (uid) => Math.max(0, metrics.height - metrics.distances.get(uid));
    for (const node of metrics.canonical) {
      if (!node.children || heightFor(node.uid) >= threshold) continue;
      const parent = metrics.parents.get(node.uid);
      if (parent === undefined || heightFor(parent) >= threshold) {
        groups.set(node.uid, policy.palette[groups.size % policy.palette.length]);
      }
    }
  }
  const base = { colour: config.branch?.colour || "#555", width: Number(config.branch?.width) || 1.5 };
  const branches = new Map();
  const stems = new Map();
  const annotations = new Map();
  const visit = (uid, inherited = {}, groupColour = null) => {
    const node = metrics.nodes.get(uid);
    const automatic = groups.get(uid) || groupColour;
    const subtree = { ...inherited, ...tree.presentation?.subtrees?.[uid]?.branch };
    const mergeStyle = { ...base, ...(automatic ? { colour: automatic } : {}), ...subtree };
    if (node.children) stems.set(uid, mergeStyle);
    const label = tree.presentation?.nodes?.[uid]?.label;
    if (label?.text) annotations.set(uid, { colour: config.textColour || "#333", fontSize: 10, fontWeight: 400, fontStyle: "normal", ...label });
    for (const child of node.children || []) {
      branches.set(child, { ...mergeStyle, ...tree.presentation?.nodes?.[child]?.branch });
      visit(child, subtree, automatic);
    }
  };
  visit(tree.root);
  const value = { base, branches, stems, annotations, policy, height: metrics.height, threshold, metric: policy.mode === "cutHeight" };
  stylesCache.set(tree, { presentation: tree.presentation, signature, value });
  return value;
}

const [markup$1, css$1] = await Promise.all([
  "<fieldset hidden>\n  <legend>Clade colouring</legend>\n  <label><span>Colour clades</span><input type=\"checkbox\" data-enabled aria-label=\"Colour clades\"></label>\n  <div data-options hidden>\n    <label><span>Cut height</span><input type=\"number\" min=\"0\" step=\"any\" data-height aria-label=\"Cut height\"><button type=\"button\" data-auto title=\"Use automatic cut height (70%)\">Auto</button></label>\n    <input type=\"range\" min=\"0\" aria-label=\"Preview cut height\">\n    <label><span>Palette</span><select aria-label=\"Clade palette\"><option value=\"custom\" hidden>Custom</option><option value=\"muted\">Muted</option><option value=\"accessible\">Colorblind friendly</option></select></label>\n    <label><span>Show cut line</span><input type=\"checkbox\" data-line aria-label=\"Show cut line\"></label>\n  </div>\n</fieldset>\n<button type=\"button\" data-reset title=\"Clear manual tree styles and labels; keep automatic clade colouring.\" disabled>Reset all manual styles</button>\n<p role=\"status\" hidden></p>\n",
  ":host { display: block; min-width: 0; color: inherit; font: inherit; }\nfieldset { min-width: 0; margin: 0 0 10px; padding: 8px; border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de)); border-radius: 6px; }\nlegend { padding: 0 4px; font-weight: 600; }\nlabel { display: flex; align-items: center; gap: 6px; margin: 6px 0; }\nlabel > span { flex: 1; min-width: 0; }\ninput[type=number] { width: 64px; min-width: 0; }\ninput[type=range] { width: 100%; }\nselect { width: 132px; min-width: 0; }\nbutton[aria-pressed=true] { font-weight: 600; }\np { color: var(--cm-danger,var(--_cm-danger, #b42318)); margin: 6px 0 0; }\n",
]);
const template$1 = createTemplate(markup$1);
const palettes = { muted: TREE_COLOURING_DEFAULTS.palette, accessible: ["#0072b2", "#e69f00", "#009e73", "#cc79a7", "#d55e00", "#56b4e9"] };
const HTMLElementBase$1 = globalThis.HTMLElement || class {};

class TreeAppearance extends HTMLElementBase$1 {
  #state;
  #height;
  #slider;
  #previewing = false;
  #cancelled = false;

  constructor() {
    super();
    const root = mountComponent(this, template$1, css$1);
    this.#height = root.querySelector("[data-height]"); this.#slider = root.querySelector('[type="range"]');
    root.querySelector("[data-enabled]").onchange = (event) => this.#change({ mode: event.target.checked ? "cutHeight" : "none" });
    this.#height.onchange = () => { if (this.#height.reportValidity()) this.#change({ cutHeight: Number(this.#height.value) }); };
    root.querySelector("[data-auto]").onclick = () => this.#change({ cutHeight: "auto" });
    root.querySelector("select").onchange = (event) => this.#change({ palette: palettes[event.target.value] || this.#state.policy.palette });
    root.querySelector("[data-line]").onchange = (event) => this.#change({ showCutLine: event.target.checked });
    root.querySelector("[data-reset]").onclick = () => this.#emit("reset-styles");
    this.#slider.oninput = () => {
      this.#cancelled = false; this.#previewing = true; this.#height.value = this.#slider.value;
      if (!this.#emit("colouring-preview", { cutHeight: Number(this.#slider.value) })) this.#cancelPreview();
    };
    this.#slider.onchange = () => { if (!this.#cancelled) this.#change({ cutHeight: Number(this.#slider.value) }); };
    this.#slider.onkeydown = (event) => { if (event.key === "Escape") { this.#cancelPreview(); event.preventDefault(); } };
    this.#slider.onpointercancel = () => this.#cancelPreview();
  }

  set state({ tree, config }) {
    if (!tree) { if (this.#previewing) this.#cancelPreview(); this.#state = null; return; }
    const { policy, height } = resolveTreePresentation(tree, config);
    let compatible = tree.distance?.kind === "mergeHeight";
    try { validateTreePresentation({ ...tree, presentation: { ...tree.presentation, colouring: { ...policy, mode: "cutHeight" } } }, config); }
    catch { compatible = false; }
    this.#state = { policy, height, compatible,
      canReset: Boolean(Object.keys(tree.presentation?.nodes || {}).length || Object.keys(tree.presentation?.subtrees || {}).length) };
    this.#previewing = false; this.#sync();
  }
  set error(message) { const status = this.shadowRoot.querySelector("p"); status.textContent = message; status.hidden = !message; }
  disconnectedCallback() { if (this.#previewing) this.#cancelPreview(); }

  #emit(type, detail) { return this.dispatchEvent(new CustomEvent(type, { detail, bubbles: true, composed: true, cancelable: true })); }
  #change(changes) { if (!this.#emit("colouring-change", changes)) this.#cancelPreview(); }
  #cancelPreview() {
    this.#cancelled = true; this.#previewing = false; this.#emit("colouring-preview", null); this.#sync();
  }
  #sync() {
    if (!this.#state) return;
    const { policy, height, compatible, canReset } = this.#state, root = this.shadowRoot;
    root.querySelector("fieldset").hidden = !compatible;
    root.querySelector("[data-enabled]").checked = policy.mode === "cutHeight";
    root.querySelector("[data-options]").hidden = policy.mode !== "cutHeight";
    const actual = policy.cutHeight === "auto" ? height * 0.7 : policy.cutHeight;
    this.#height.value = actual;
    Object.assign(this.#slider, { max: String(Math.max(height, actual) || 1), step: String(height / 1000 || 0.001), value: actual });
    root.querySelector("[data-auto]").setAttribute("aria-pressed", String(policy.cutHeight === "auto"));
    const palette = Object.keys(palettes).find((key) => JSON.stringify(palettes[key]) === JSON.stringify(policy.palette)) || "custom";
    root.querySelector('[value="custom"]').hidden = palette !== "custom";
    root.querySelector("select").value = palette;
    root.querySelector("[data-line]").checked = policy.showCutLine;
    root.querySelector("[data-reset]").disabled = !canReset;
  }
}

function defineTreeAppearance() {
  if (!customElements.get("clinker-tree-appearance")) customElements.define("clinker-tree-appearance", TreeAppearance);
}

uiAppearanceStyles(":where(.cm-chart,.cm-ui)") + lightControlStyles +
  ".cm-chart { background:var(--cm-background,var(--_cm-background)); }\n.cm-chart > :is(svg.clusterMap,canvas.clusterMapCanvas,canvas.clusterMapWebGpuOverlay) { z-index:1; background:transparent; color-scheme:light; }\n.cm-chart > .clusterMapFigureSurface { z-index:0; box-shadow:0 2px 12px rgba(0,0,0,.2); outline:1px solid var(--cm-border,var(--_cm-border)); }\n.cm-chart > .clusterMapFigureDimensions { z-index:2; color:var(--cm-muted,var(--_cm-muted)); font:var(--cm-font,var(--_cm-font)); line-height:20px; white-space:nowrap; }\n.cm-chart > canvas.clusterMapMinimap { color-scheme:light; }\n\n.cm-chart > :is(svg.clusterMap,canvas.clusterMapCanvas) { position:relative; }\n.cm-chart > :is(canvas.clusterMapCanvas,canvas.clusterMapWebGpuOverlay) { display:block; width:100%; height:100%; }\n.cm-chart > canvas.clusterMapCanvas { outline:none; cursor:grab; }\n.cm-chart > canvas.clusterMapWebGpuOverlay { position:absolute; inset:0; pointer-events:none; }\n.cm-chart > canvas.clusterMapMinimap { position:absolute; z-index:2; display:block; box-sizing:border-box; background:white; box-shadow:0 1px 4px rgba(0,0,0,.25); cursor:grab; touch-action:none; }\n.cm-chart > :is(.clusterMapFigureSurface,.clusterMapFigureDimensions) { position:absolute; left:0; top:0; pointer-events:none; }\n.cm-chart > .clusterMapFigureSurface { transform-origin:0 0; }\n.cm-chart > .clusterMapFigureSurface--frameless { outline:none; box-shadow:none; }\n.cm-chart > input.colourPicker { position:absolute; opacity:0; }\n.cm-chart > div.tooltip { position:absolute; pointer-events:none; z-index:4; box-sizing:border-box; padding:8px; background:var(--cm-surface,var(--_cm-surface)); border:1px solid var(--cm-border,var(--_cm-border)); border-radius:4px; box-shadow:0 2px 8px rgba(0,0,0,.2); }\n";

/** Figure colours are independent of the interface/system theme. */
function figureConfig(config) {
  if (config.figureStyle) return config;
  const theme = { background: config.plot?.background ?? "#ffffff", foreground: "#000000",
    outline: "#000000", branch: "#555555", guide: "#999999",
    muted: "#687586", missing: "#eeeeee" };
  return { ...config, figureStyle: theme,
    tree: { ...config.tree, textColour: theme.foreground,
      branch: { ...config.tree?.branch, colour: config.tree?.branch?.colour ?? theme.branch },
      leafGuide: { ...config.tree?.leafGuide, colour: config.tree?.leafGuide?.colour ?? theme.guide } },
    gene: { ...config.gene, shape: { ...config.gene?.shape, stroke: config.gene?.shape?.stroke ?? theme.outline } },
    locus: { ...config.locus, trackBar: { ...config.locus?.trackBar, colour: config.locus?.trackBar?.colour ?? theme.outline } },
    scaleBar: { ...config.scaleBar, colour: config.scaleBar?.colour ?? theme.foreground },
  };
}

const [markup, css] = await Promise.all([
  "<div class=\"toolbar\" part=\"toolbar\"><clinker-search-field data-appearance-filter label=\"Find an appearance setting\" placeholder=\"Find a setting\"></clinker-search-field></div>\n<div class=\"settings\" part=\"settings\">\n  <fieldset data-appearance-section=\"plot\">\n    <legend>Plot</legend>\n    <label><span>Show figure frame</span><input type=\"checkbox\" data-config=\"plot.showFigureFrame\" aria-label=\"Show figure border and dimensions\"></label>\n    <label><span>Figure background</span><input type=\"color\" data-config=\"plot.background\" aria-label=\"Figure background\"></label>\n    <label><span>Enable annotation panels</span><input type=\"checkbox\" data-config=\"annotations.enabled\" aria-label=\"Enable annotation panels\"></label>\n  </fieldset>\n  <fieldset data-appearance-section=\"legend\">\n    <legend>Legend</legend>\n    <label><span>Show legend</span><input type=\"checkbox\" data-config=\"legend.show\" aria-label=\"Show legend\"></label>\n    <label><span>Position</span><select data-config=\"legend.position\" aria-label=\"Position\"><option value=\"right\">right</option><option value=\"bottom\">bottom</option></select></label>\n    <label><span>Columns</span><input type=\"number\" step=\"any\" data-config=\"legend.columns\" aria-label=\"Columns\"></label>\n    <details data-appearance-section=\"legend\">\n      <summary>Advanced</summary>\n      <label><span>Column width</span><input type=\"number\" step=\"any\" data-config=\"legend.columnWidth\" aria-label=\"Column width\"></label>\n      <label><span>Entry height</span><input type=\"number\" step=\"any\" data-config=\"legend.entryHeight\" aria-label=\"Entry height\"></label>\n      <label><span>Font size</span><input type=\"number\" step=\"any\" data-config=\"legend.fontSize\" aria-label=\"Font size\"></label>\n      <label><span>Subtitle size</span><input type=\"number\" step=\"any\" data-config=\"legend.subtitleFontSize\" aria-label=\"Subtitle size\"></label>\n      <label><span>Left margin</span><input type=\"number\" step=\"any\" data-config=\"legend.marginLeft\" aria-label=\"Left margin\"></label>\n      <label><span>Top margin</span><input type=\"number\" step=\"any\" data-config=\"legend.marginTop\" aria-label=\"Top margin\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"links\">\n    <legend>Links</legend>\n    <label><span>Show links</span><input type=\"checkbox\" data-config=\"link.show\" aria-label=\"Show links\"></label>\n    <label><span>Draw as lines</span><input type=\"checkbox\" data-config=\"link.asLine\" aria-label=\"Draw as lines\"></label>\n    <label><span>Best match only</span><input type=\"checkbox\" data-config=\"link.bestOnly\" aria-label=\"Best match only\"></label>\n    <label><span>Minimum identity</span><input type=\"number\" step=\"any\" data-config=\"link.threshold\" aria-label=\"Minimum identity\"></label>\n    <details data-appearance-section=\"links\">\n      <summary>Advanced</summary>\n      <label><span>Straight links</span><input type=\"checkbox\" data-config=\"link.straight\" aria-label=\"Straight links\"></label>\n      <label><span>Use group colour</span><input type=\"checkbox\" data-config=\"link.groupColour\" aria-label=\"Use group colour\"></label>\n      <label><span>Stroke width</span><input type=\"number\" step=\"any\" data-config=\"link.strokeWidth\" aria-label=\"Stroke width\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"link-labels\">\n    <legend>Link labels</legend>\n    <label><span>Show link labels</span><input type=\"checkbox\" data-config=\"link.label.show\" aria-label=\"Show link labels\"></label>\n    <details data-appearance-section=\"link-labels\">\n      <summary>Advanced</summary>\n      <label><span>Label background</span><input type=\"checkbox\" data-config=\"link.label.background\" aria-label=\"Label background\"></label>\n      <label><span>Font size</span><input type=\"number\" step=\"any\" data-config=\"link.label.fontSize\" aria-label=\"Font size\"></label>\n      <label><span>Position</span><input type=\"number\" step=\"any\" data-config=\"link.label.position\" aria-label=\"Position\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"tree\">\n    <legend>Tree</legend>\n    <label><span>Show tree</span><input type=\"checkbox\" data-config=\"tree.show\" aria-label=\"Show tree\"></label>\n    <label><span>Show tree nodes</span><input type=\"checkbox\" data-config=\"tree.node.show\" aria-label=\"Show tree nodes\"></label>\n    <label><span>Tree width</span><input type=\"number\" step=\"any\" data-config=\"tree.width\" aria-label=\"Tree width\"></label>\n    <label><span>Use branch lengths</span><input type=\"checkbox\" data-config=\"tree.useBranchLengths\" aria-label=\"Use branch lengths\"></label>\n    <label><span>Branch stroke width</span><input type=\"number\" step=\"any\" data-config=\"tree.branch.width\" aria-label=\"Branch stroke width\"></label>\n    <label><span>Branch colour</span><input type=\"color\" data-config=\"tree.branch.colour\" aria-label=\"Branch colour\"></label>\n    <label><span>Show support values</span><input type=\"checkbox\" data-config=\"tree.support.show\" aria-label=\"Show support values\"></label>\n    <label><span>Show tree scale bar</span><input type=\"checkbox\" data-config=\"tree.scaleBar.show\" aria-label=\"Show tree scale bar\"></label>\n    <label><span>Show dashed leaf guides</span><input type=\"checkbox\" data-config=\"tree.leafGuide.show\" aria-label=\"Show dashed leaf guides\"></label>\n    <details data-appearance-section=\"tree\">\n      <summary>Advanced</summary>\n      <label><span>Support font size</span><input type=\"number\" step=\"any\" data-config=\"tree.support.fontSize\" aria-label=\"Support font size\"></label>\n      <label><span>Support precision</span><input type=\"number\" step=\"any\" data-config=\"tree.support.precision\" aria-label=\"Support precision\"></label>\n      <label><span>Leaf guide colour</span><input type=\"color\" data-config=\"tree.leafGuide.colour\" aria-label=\"Leaf guide colour\"></label>\n      <label><span>Leaf guide width</span><input type=\"number\" step=\"any\" data-config=\"tree.leafGuide.width\" aria-label=\"Leaf guide width\"></label>\n      <label><span>Leaf guide dash pattern</span><input type=\"text\" data-config=\"tree.leafGuide.dash\" aria-label=\"Leaf guide dash pattern\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"clusters-loci\">\n    <legend>Clusters &amp; loci</legend>\n    <label><span>Align cluster labels</span><input type=\"checkbox\" data-config=\"cluster.alignLabels\" aria-label=\"Align cluster labels\"></label>\n    <label><span>Hide locus coordinates</span><input type=\"checkbox\" data-config=\"cluster.hideLocusCoordinates\" aria-label=\"Hide locus coordinates\"></label>\n    <label><span>Show trim indicators</span><input type=\"checkbox\" data-config=\"locus.trimIndicators.show\" aria-label=\"Show trim indicators\"></label>\n    <label><span>Cluster spacing</span><input type=\"number\" step=\"any\" data-config=\"cluster.spacing\" aria-label=\"Cluster spacing\"></label>\n    <label><span>Locus spacing</span><input type=\"number\" step=\"any\" data-config=\"locus.spacing\" aria-label=\"Locus spacing\"></label>\n    <details data-appearance-section=\"clusters-loci\">\n      <summary>Advanced</summary>\n      <label><span>Cluster name size</span><input type=\"number\" step=\"any\" data-config=\"cluster.nameFontSize\" aria-label=\"Cluster name size\"></label>\n      <label><span>Locus label size</span><input type=\"number\" step=\"any\" data-config=\"cluster.lociFontSize\" aria-label=\"Locus label size\"></label>\n      <label><span>Trim indicator size</span><input type=\"number\" step=\"any\" data-config=\"locus.trimIndicators.fontSize\" aria-label=\"Trim indicator size\"></label>\n      <label><span>Trim indicator weight</span><input type=\"number\" step=\"any\" data-config=\"locus.trimIndicators.fontWeight\" aria-label=\"Trim indicator weight\"></label>\n      <label><span>Trim indicator border</span><input type=\"number\" step=\"any\" data-config=\"locus.trimIndicators.borderWidth\" aria-label=\"Trim indicator border\"></label>\n      <label><span>Trim indicator border opacity</span><input type=\"number\" step=\"any\" data-config=\"locus.trimIndicators.borderOpacity\" aria-label=\"Trim indicator border opacity\"></label>\n      <label><span>Track colour</span><input type=\"color\" data-config=\"locus.trackBar.colour\" aria-label=\"Track colour\"></label>\n      <label><span>Track stroke</span><input type=\"number\" step=\"any\" data-config=\"locus.trackBar.stroke\" aria-label=\"Track stroke\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"genes\">\n    <legend>Genes</legend>\n    <label><span>Show gene labels</span><input type=\"checkbox\" data-config=\"gene.label.show\" aria-label=\"Show gene labels\"></label>\n    <label><span>Label field</span><select data-config=\"gene.label.name\" aria-label=\"Label field\"><option value=\"uid\">uid</option><option value=\"name\">name</option><option value=\"label\">label</option></select></label>\n    <label><span>Label position</span><select data-config=\"gene.label.position\" aria-label=\"Label position\"><option value=\"top\">top</option><option value=\"bottom\">bottom</option></select></label>\n    <details data-appearance-section=\"genes\">\n      <summary>Advanced</summary>\n      <label><span>Label anchor</span><select data-config=\"gene.label.anchor\" aria-label=\"Label anchor\"><option value=\"start\">start</option><option value=\"middle\">middle</option><option value=\"end\">end</option></select></label>\n      <label><span>Label size</span><input type=\"number\" step=\"any\" data-config=\"gene.label.fontSize\" aria-label=\"Label size\"></label>\n      <label><span>Label rotation</span><input type=\"number\" step=\"any\" data-config=\"gene.label.rotation\" aria-label=\"Label rotation\"></label>\n      <label><span>Label spacing</span><input type=\"number\" step=\"any\" data-config=\"gene.label.spacing\" aria-label=\"Label spacing\"></label>\n      <label><span>Label start</span><input type=\"number\" step=\"any\" data-config=\"gene.label.start\" aria-label=\"Label start\"></label>\n      <label><span>Body height</span><input type=\"number\" step=\"any\" data-config=\"gene.shape.bodyHeight\" aria-label=\"Body height\"></label>\n      <label><span>Tip height</span><input type=\"number\" step=\"any\" data-config=\"gene.shape.tipHeight\" aria-label=\"Tip height\"></label>\n      <label><span>Tip length</span><input type=\"number\" step=\"any\" data-config=\"gene.shape.tipLength\" aria-label=\"Tip length\"></label>\n      <label><span>Stroke colour</span><input type=\"color\" data-config=\"gene.shape.stroke\" aria-label=\"Stroke colour\"></label>\n      <label><span>Stroke width</span><input type=\"number\" step=\"any\" data-config=\"gene.shape.strokeWidth\" aria-label=\"Stroke width\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"scale-bars\">\n    <legend>Scale bars</legend>\n    <label><span>Show scale bar</span><input type=\"checkbox\" data-config=\"scaleBar.show\" aria-label=\"Show scale bar\"></label>\n    <label><span>Base pairs</span><input type=\"number\" step=\"any\" data-config=\"scaleBar.basePair\" aria-label=\"Base pairs\"></label>\n    <label><span>Show identity bar</span><input type=\"checkbox\" data-config=\"colourBar.show\" aria-label=\"Show identity bar\"></label>\n    <label><span>Minimum source</span><select data-config=\"colourBar.domain.minMode\" aria-label=\"Minimum source\"><option value=\"fixed\">Custom value</option><option value=\"data\">Lowest identity in data</option></select></label>\n    <label><span>Minimum identity</span><input type=\"number\" step=\"any\" data-config=\"colourBar.domain.min\" aria-label=\"Minimum identity\"></label>\n    <label><span>Maximum source</span><select data-config=\"colourBar.domain.maxMode\" aria-label=\"Maximum source\"><option value=\"fixed\">Custom value</option><option value=\"data\">Highest identity in data</option></select></label>\n    <label><span>Maximum identity</span><input type=\"number\" step=\"any\" data-config=\"colourBar.domain.max\" aria-label=\"Maximum identity\"></label>\n    <details data-appearance-section=\"scale-bars\">\n      <summary>Advanced</summary>\n      <label><span>Height</span><input type=\"number\" step=\"any\" data-config=\"scaleBar.height\" aria-label=\"Height\"></label>\n      <label><span>Font size</span><input type=\"number\" step=\"any\" data-config=\"scaleBar.fontSize\" aria-label=\"Font size\"></label>\n      <label><span>Stroke</span><input type=\"number\" step=\"any\" data-config=\"scaleBar.stroke\" aria-label=\"Stroke\"></label>\n      <label><span>Colour</span><input type=\"color\" data-config=\"scaleBar.colour\" aria-label=\"Colour\"></label>\n      <label><span>Scale margin</span><input type=\"number\" step=\"any\" data-config=\"scaleBar.marginTop\" aria-label=\"Scale margin\"></label>\n      <label><span>Identity width</span><input type=\"number\" step=\"any\" data-config=\"colourBar.width\" aria-label=\"Identity width\"></label>\n      <label><span>Identity height</span><input type=\"number\" step=\"any\" data-config=\"colourBar.height\" aria-label=\"Identity height\"></label>\n      <label><span>Identity font size</span><input type=\"number\" step=\"any\" data-config=\"colourBar.fontSize\" aria-label=\"Identity font size\"></label>\n      <label><span>Identity margin</span><input type=\"number\" step=\"any\" data-config=\"colourBar.marginTop\" aria-label=\"Identity margin\"></label>\n    </details>\n  </fieldset>\n  <fieldset data-appearance-section=\"plot-navigation\">\n    <legend>Plot &amp; navigation</legend>\n    <label><span>Renderer</span><select data-config=\"plot.renderer\" aria-label=\"Renderer\"><option value=\"svg\">svg</option><option value=\"canvas\">canvas</option><option value=\"webgpu\">webgpu</option></select></label>\n    <label><span>Show minimap</span><input type=\"checkbox\" data-config=\"plot.minimap.show\" aria-label=\"Show minimap\"></label>\n    <label><span>Minimum zoom</span><input type=\"number\" step=\"any\" data-config=\"plot.minZoom\" aria-label=\"Minimum zoom\"></label>\n    <label><span>Maximum zoom</span><input type=\"number\" step=\"any\" data-config=\"plot.maxZoom\" aria-label=\"Maximum zoom\"></label>\n    <details data-appearance-section=\"plot-navigation\">\n      <summary>Advanced</summary>\n      <label><span>Transition (ms)</span><input type=\"number\" step=\"any\" data-config=\"plot.transitionDuration\" aria-label=\"Transition (ms)\"></label>\n      <label><span>Scale factor</span><input type=\"number\" step=\"any\" data-config=\"plot.scaleFactor\" aria-label=\"Scale factor\"></label>\n      <label><span>Scale genes</span><input type=\"checkbox\" data-config=\"plot.scaleGenes\" aria-label=\"Scale genes\"></label>\n      <label><span>Font family</span><input type=\"text\" data-config=\"plot.fontFamily\" aria-label=\"Font family\"></label>\n      <label><span>Minimap links</span><input type=\"checkbox\" data-config=\"plot.minimap.showLinks\" aria-label=\"Minimap links\"></label>\n      <label><span>Minimap width</span><input type=\"number\" step=\"any\" data-config=\"plot.minimap.width\" aria-label=\"Minimap width\"></label>\n      <label><span>Minimap height</span><input type=\"number\" step=\"any\" data-config=\"plot.minimap.height\" aria-label=\"Minimap height\"></label>\n      <label><span>Minimap margin</span><input type=\"number\" step=\"any\" data-config=\"plot.minimap.margin\" aria-label=\"Minimap margin\"></label>\n    </details>\n  </fieldset>\n  <p class=\"empty\" hidden></p>\n  <clinker-tree-appearance hidden></clinker-tree-appearance>\n  <p class=\"error\" role=\"alert\" hidden></p>\n</div>\n",
  ":host { display: flex; flex: 1; flex-direction: column; min-width: 0; min-height: 0; gap: 8px; }\n.toolbar { display: flex; flex-shrink: 0; }\n.settings { display: grid; grid-template-columns: repeat(auto-fit, minmax(min(100%, 210px), 1fr)); gap: 10px 20px; align-content: start; overflow: auto; padding: 4px 0; min-height: 0; }\nfieldset { min-width: 0; margin: 0; padding: 8px; border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de)); border-radius: 6px; }\nlegend { padding: 0 4px; font-weight: 600; }\nlabel { display: flex; align-items: center; gap: 6px; margin: 6px 0; }\nlabel > span { min-width: 0; flex: 1; }\nselect, input:not([type=checkbox]):not([type=color]) { width: 88px; min-width: 0; }\ndetails { margin-top: 8px; border-top: 1px solid var(--cm-border,var(--_cm-border, #e5e9ef)); }\nsummary { padding: 7px 0 1px; color: var(--cm-muted,var(--_cm-muted, #687586)); cursor: pointer; font-size: 12px; font-weight: 600; }\nsummary:hover { color: var(--cm-accent,var(--_cm-accent, #1677ff)); }\n.empty { grid-column: 1 / -1; padding: 16px; color: var(--cm-muted,var(--_cm-muted, #667386)); text-align: center; }\n.error { grid-column: 1 / -1; color: var(--cm-danger,var(--_cm-danger, #b42318)); margin: 0; }\n",
]);
const template = createTemplate(markup);
const HTMLElementBase = globalThis.HTMLElement || class {};

class ClinkerAppearancePanelElement extends HTMLElementBase {
  static observedAttributes = ["hidden"];
  #chart = null;
  #unsubscribe = null;
  #config;
  #tree;
  #runEdit;
  #hasTree = false;
  #filter = "";
  #expanded = new Set();
  #sections;
  #navigation = 0;

  constructor() {
    super();
    const root = mountComponent(this, template, editorThemeStyles + css);
    this.#tree = root.querySelector("clinker-tree-appearance");
    this.#runEdit = createActionRunner((message) => {
      const status = root.querySelector(".error"); status.textContent = message; status.hidden = !message;
    });
    const runTreeEdit = createActionRunner((message) => { this.#tree.error = message; });
    const runPreview = createActionRunner((message) => { if (message) this.#tree.error = message; });
    const colouring = (changes) => [{ type: "tree.presentation.update", scope: "tree", changes: { colouring: changes } }];
    this.#tree.addEventListener("colouring-change", (event) => { if (!runTreeEdit(() => this.#chart.patch(colouring(event.detail)))) event.preventDefault(); });
    this.#tree.addEventListener("colouring-preview", (event) => { if (!runPreview(() => this.#chart?.previewTreePresentation(event.detail ? colouring(event.detail) : null))) event.preventDefault(); });
    this.#tree.addEventListener("reset-styles", (event) => { if (!runTreeEdit(() => this.#chart.patch([{ type: "tree.presentation.reset", scope: "manual" }]))) event.preventDefault(); });
    this.#sections = [...root.querySelectorAll("fieldset")].map((box) => ({ box, key: box.dataset.appearanceSection, details: box.querySelector("details"),
      fields: [...box.querySelectorAll("[data-config]")].map((control) => ({ control, row: control.closest("label"), advanced: Boolean(control.closest("details")),
        search: `${box.querySelector("legend").textContent} ${control.getAttribute("aria-label")} ${control.dataset.config}`.toLowerCase() })) }));
    const search = root.querySelector("clinker-search-field");
    search.addEventListener("input", () => { this.#filter = search.value.trim().toLowerCase(); this.#sync(); });
    root.addEventListener("click", (event) => {
      const details = event.target.closest("summary")?.parentElement;
      if (!details || event.defaultPrevented || this.#filter) return;
      const key = details.dataset.appearanceSection;
      if (details.open) this.#expanded.delete(key); else this.#expanded.add(key);
    });
    root.addEventListener("change", (event) => {
      const input = event.target;
      if (!input.matches("[data-config]")) return;
      if (!input.reportValidity()) return;
      const value = input.type === "checkbox" ? input.checked : input.type === "number" ? Number(input.value) : input.value;
      const config = input.dataset.config.split(".").reduceRight((value, key) => ({ [key]: value }), value);
      if (!this.#runEdit(() => this.#chart.config(config))) this.#syncField(input, true);
    });
    upgradeProperty(this, "chart");
  }

  get chart() { return this.#chart; }
  get sections() { return this.#sections.map(({ key }) => key); }

  async openSection(section = null, { focus = true } = {}) {
    const target = section === null ? this.#sections[0] : this.#sections.find(({ key }) => key === section);
    if (!target || !this.isConnected || !this.#chart) return false;
    const navigation = ++this.#navigation;
    this.hidden = false;
    this.#filter = "";
    this.shadowRoot.querySelector("clinker-search-field").value = "";
    if (target.details) this.#expanded.add(target.details.dataset.appearanceSection);
    this.#render();
    await new Promise((resolve) => this.ownerDocument.defaultView.requestAnimationFrame(resolve));
    // A later navigation, close or tab switch must not steal focus back.
    if (navigation !== this.#navigation || !this.isConnected || !target.box.getClientRects().length) return false;
    const control = target.fields.find(({ control }) => !control.disabled)?.control;
    (control || target.box).scrollIntoView({ block: "nearest" });
    if (focus) control?.focus({ preventScroll: true });
    return true;
  }
  set chart(value) {
    if (value === this.#chart) return;
    this.#navigation++;
    this.#chart?.previewTreePresentation?.(null);
    this.#unsubscribe?.(); this.#unsubscribe = null;
    this.#chart = value || null; this.#subscribe(); this.#render();
  }
  connectedCallback() {
    this.classList.add("cm-ui", "cm-data-ui");
    this.#subscribe(); this.#render();
  }
  disconnectedCallback() {
    this.#navigation++;
    this.#chart?.previewTreePresentation?.(null);
    this.#unsubscribe?.(); this.#unsubscribe = null;
  }
  attributeChangedCallback() { if (this.isConnected && !this.hidden) this.#render(); }
  #subscribe() {
    if (this.#unsubscribe || !this.isConnected || !this.#chart?.on) return;
    this.#unsubscribe = this.#chart.on("change", (change) => {
      if (["config.change", "theme.change", "data.replace", "state.replace"].includes(change.type) ||
          (change.type === "data.apply" && change.operations?.some((op) => op.type.startsWith("tree.")))) this.#render();
    });
  }
  #render() {
    if (!this.isConnected) return;
    const config = this.#chart?.config?.();
    applyEditorTheme(this, this.#chart);
    this.shadowRoot.querySelector(".toolbar").hidden = !config;
    this.shadowRoot.querySelector(".settings").hidden = !config;
    if (!config) { this.#config = null; this.#hasTree = false; this.#tree.state = { tree: null }; return; }
    if (this.hidden) return;
    this.#config = figureConfig(config);
    const tree = this.#chart.data()?.tree;
    this.#hasTree = Boolean(tree); this.#tree.state = { tree, config: this.#config.tree };
    this.#sync();
  }

  #syncField(control, force = false) {
    const path = control.dataset.config, value = path.split(".").reduce((current, key) => current?.[key], this.#config);
    if (control.type === "checkbox") control.checked = Boolean(value);
    else if (force || control !== this.shadowRoot.activeElement) control.value = control.type === "color" ? colourHex(value) : value ?? "";
    const bound = path === "colourBar.domain.min" ? "min" : path === "colourBar.domain.max" ? "max" : null;
    control.disabled = Boolean(bound && this.#config?.colourBar?.domain?.[`${bound}Mode`] === "data");
  }

  #sync() {
    let visible = 0;
    for (const { box, details, fields } of this.#sections) {
      for (const { control, row, search } of fields) {
        this.#syncField(control); row.hidden = Boolean(this.#filter && !search.includes(this.#filter));
      }
      box.hidden = fields.every(({ row }) => row.hidden);
      if (!box.hidden) visible++;
      if (details) {
        details.hidden = !fields.some(({ advanced, row }) => advanced && !row.hidden);
        details.open = Boolean(this.#filter) || this.#expanded.has(details.dataset.appearanceSection);
      }
    }
    const empty = this.shadowRoot.querySelector(".empty");
    empty.hidden = visible > 0; empty.textContent = `No appearance settings match “${this.#filter}”.`;
    this.#tree.hidden = !this.#hasTree || Boolean(this.#filter && !/tree|clade|colour|color|annotation|branch|cut|palette/.test(this.#filter));
  }
}

function defineAppearancePanel() {
  defineSearchField(); defineTreeAppearance();
  if (!customElements.get("clinker-appearance-panel")) customElements.define("clinker-appearance-panel", ClinkerAppearancePanelElement);
}

export { defineAppearancePanel as d };
