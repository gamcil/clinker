import { d as defineActionMenu, i as installDataEditorStyle, a as installStyle } from './editor/annotation-editor-3c63735d.mjs';
import { c as createTemplate, m as mountComponent, a as controlStyles, u as upgradeProperty } from './editor/search-field-595834d9.mjs';
import { a as applyEditorTheme, e as editorThemeStyles } from './editor/theme-8bd47f2f.mjs';
import { d as defineAppearancePanel } from './editor/appearance-panel-b4ed47ff.mjs';
import { b as bindOptionButtons, d as defineDataPanel } from './editor/data-panel-a0484e62.mjs';

const editable = "input, textarea, select, [contenteditable]:not([contenteditable='false'])";

function eventMatches(event, selector) {
  const path = event.composedPath?.() || [event.target];
  return path.some((node) => node?.closest?.(selector) || node?.isContentEditable);
}

const isEditingEvent = (event) => eventMatches(event, editable);

const toolbarTemplate = createTemplate("<strong></strong>\n<button type=\"button\" data-action=\"focus\">Focus</button>\n<button type=\"button\" data-action=\"flip\">Flip</button>\n<button type=\"button\" data-action=\"trim\">Trim</button>\n<clinker-action-menu label=\"Loci…\" data-loci-menu hidden>\n  <button type=\"button\" data-action=\"restore\">Restore full loci</button>\n</clinker-action-menu>\n<clinker-action-menu label=\"More…\" data-more-menu>\n  <button type=\"button\" data-action=\"edit\">Edit data…</button>\n  <button type=\"button\" data-action=\"remove\"></button>\n</clinker-action-menu>\n<button type=\"button\" data-action=\"clear\">Clear</button>\n");
const css$2 = ".cm-plot-selection-toolbar { position: absolute; z-index: 2; top: 12px; left: 12px; display: flex; align-items: center; gap: 6px; padding: 7px 8px; border: 1px solid var(--cm-border,var(--_cm-border,#c7d9f1)); border-radius: 8px; background: var(--cm-surface,var(--_cm-surface,rgba(255,255,255,.94))); box-shadow: 0 3px 12px rgba(23,36,54,.12); color: var(--cm-foreground,var(--_cm-foreground,#273444)); font: var(--cm-font,var(--_cm-font)); }\n.cm-plot-selection-toolbar[hidden] { display: none; }\n.cm-plot-selection-toolbar strong { padding: 4px 6px; border-radius: 5px; background: var(--cm-hover,var(--_cm-hover,#e7f0ff)); color: var(--cm-accent,var(--_cm-accent,#184a86)); }\n";

let activeToolbar = null;

/** Mount contextual plot actions without requiring the full data editor. */
function mountPlotSelectionToolbar(chart, container, { editor = null } = {}) {
  if (!chart?.selection || !chart?.flipLoci || !chart?.trimLoci || !chart?.restoreLoci) throw new TypeError("A mounted ClusterMap chart is required.");
  if (!(container instanceof Element)) throw new TypeError("A plot container element is required.");
  const document = container.ownerDocument;
  const scope = container.getRootNode();
  defineActionMenu(); installDataEditorStyle(scope); installStyle(scope, "clinker-plot-selection-style", css$2);
  const root = document.createElement("div");
  root.className = "cm-ui cm-data-ui cm-plot-selection-toolbar";
  root.hidden = true;
  root.append(toolbarTemplate(document));
  const button = (action) => root.querySelector(`[data-action="${action}"]`);
  const flip = button("flip"), trim = button("trim"), restore = button("restore");
  const lociMenu = root.querySelector("[data-loci-menu]"), moreMenu = root.querySelector("[data-more-menu]");
  let clusterActions;
  const lociFor = ({ kind, ids }) => kind === "loci" ? ids :
    (chart.visibleData()?.clusters || []).filter((cluster) => ids.includes(cluster.uid)).flatMap((cluster) => cluster.loci.map((locus) => locus.uid));
  const sync = () => {
    applyEditorTheme(root, chart);
    const { kind, ids } = chart.selection(), clusters = kind === "clusters";
    root.hidden = !ids.length;
    const noun = ids.length === 1 ? (clusters ? "cluster" : "locus") : kind;
    root.querySelector("strong").textContent = `${ids.length} ${noun} selected`;
    flip.textContent = clusters ? "Flip each locus" : "Flip"; trim.textContent = clusters ? "Trim each locus" : "Trim";
    button("edit").hidden = !editor; button("remove").textContent = `Remove selected ${kind}`;
    if (clusterActions !== clusters) {
      if (clusters) lociMenu.prepend(flip, trim, restore);
      else { root.insertBefore(flip, lociMenu); root.insertBefore(trim, lociMenu); moreMenu.prepend(restore); }
      lociMenu.open = moreMenu.open = false;
      lociMenu.hidden = !clusters; clusterActions = clusters;
    }
    if (root.hidden) lociMenu.open = moreMenu.open = false;
  };
  const clear = () => chart.selection({ kind: chart.selection().kind, ids: [] });
  const activate = () => { activeToolbar = root; };
  root.addEventListener("click", (event) => {
    const action = event.target.closest("[data-action]")?.dataset.action;
    if (!action) return;
    const selection = chart.selection(), loci = lociFor(selection);
    if (action === "flip") chart.flipLoci(loci);
    else if (action === "trim") chart.trimLoci(loci);
    else if (action === "restore") chart.restoreLoci(loci);
    else if (action === "clear") clear();
    else if (action === "edit") editor?.editSelection(selection);
    else if (action === "focus") chart.focus({ [selection.kind]: selection.ids });
    else if (action === "remove" && document.defaultView.confirm(`Remove ${selection.ids.length} ${selection.kind} and their contents from the plot? This can be undone.`)) {
      chart.remove({ [selection.kind]: selection.ids });
    }
    for (const menu of root.querySelectorAll("clinker-action-menu[open]")) menu.open = false;
  });
  const onKeyDown = (event) => {
    if (activeToolbar !== root || event.key !== "Escape" || event.defaultPrevented || !chart.selection().ids.length) return;
    if (isEditingEvent(event)) return;
    event.preventDefault(); clear();
  };
  container.addEventListener("pointerdown", activate);
  container.addEventListener("focusin", activate);
  document.addEventListener("keydown", onKeyDown);
  const unsubscribe = chart.on("change", (change) => {
    if (["theme.change", "config.change"].includes(change.type)) applyEditorTheme(root, chart);
    if (["selection.change", "state.replace", "data.apply", "data.replace"].includes(change.type)) sync();
  });
  container.append(root);
  sync();
  return () => {
    unsubscribe();
    if (activeToolbar === root) activeToolbar = null;
    container.removeEventListener("pointerdown", activate);
    container.removeEventListener("focusin", activate);
    document.removeEventListener("keydown", onKeyDown);
    root.remove();
  };
}

const [markup$1, css$1] = await Promise.all([
  "<div part=\"tab-list\" role=\"tablist\">\n  <template><button type=\"button\" role=\"tab\"></button></template>\n</div>\n<slot name=\"status\"></slot>\n<div part=\"panels\">\n  <template><div role=\"tabpanel\" tabindex=\"0\"><slot></slot></div></template>\n</div>\n",
  ":host { display: flex; flex: 1; flex-direction: column; min-height: 0; color: inherit; font: inherit; }\n:host([hidden]), [hidden] { display: none !important; }\n[role=tablist] { display: flex; flex-shrink: 0; gap: 4px; margin-bottom: 10px; border-bottom: 1px solid var(--cm-border,var(--_cm-border, #e1e5ea)); }\nbutton {\n  min-height: var(--cm-control-height,var(--_cm-control-height,30px));\n  border: 0;\n  border-bottom: 2px solid transparent;\n  border-radius: 0;\n  background: transparent;\n}\nbutton[aria-selected=true] { border-bottom-color: var(--cm-accent,var(--_cm-accent, #1677ff)); font-weight: 600; }\nbutton:hover:not(:disabled) { background: var(--cm-hover,var(--_cm-hover, #e3ecf8)); }\n[part=panels], [role=tabpanel] { display: flex; flex: 1; flex-direction: column; min-height: 0; }\n[role=tabpanel]:focus-visible { outline: 2px solid var(--cm-focus-ring,var(--_cm-focus-ring, #1677ff)); outline-offset: -2px; }\nslot { display: contents; }\n",
]);
const template = createTemplate(markup$1);
const HTMLElementBase$1 = globalThis.HTMLElement || class {};
let nextId = 0;

class Tabs extends HTMLElementBase$1 {
  static observedAttributes = ["value", "label"];
  #prefix = `clinker-tabs-${++nextId}`;
  #list;
  #content;
  #buttonTemplate;
  #panelTemplate;
  #tabs = [];
  #observer;

  constructor() {
    super();
    const root = mountComponent(this, template, css$1);
    this.#list = root.querySelector("[role=tablist]"); this.#content = root.querySelector("[part=panels]");
    this.#buttonTemplate = this.#list.querySelector("template"); this.#panelTemplate = this.#content.querySelector("template");
    bindOptionButtons(this.#list, (button) => {
      if (!button || button.disabled || button.value === this.value) return;
      this.value = button.value;
      this.dispatchEvent(new CustomEvent("value-change", { detail: { value: this.value }, bubbles: true, composed: true }));
    });
    this.#observer = new this.ownerDocument.defaultView.MutationObserver((records) => {
      if (records.some((record) => record.type === "childList" ? record.target === this : record.target.parentElement === this)) this.#mountPanels();
    });
  }

  get value() { return this.getAttribute("value") || ""; }
  set value(value) { if (this.value !== String(value)) this.setAttribute("value", String(value)); }
  attributeChangedCallback() { this.#sync(); }
  connectedCallback() {
    this.#mountPanels();
    this.#observer.observe(this, { childList: true, subtree: true, attributes: true, attributeFilter: ["name", "label", "disabled"] });
  }
  disconnectedCallback() { this.#observer.disconnect(); }

  #mountPanels() {
    this.#tabs = [...this.children].filter((panel) => panel.hasAttribute("name")).map((panel, index) => {
      const button = this.ownerDocument.importNode(this.#buttonTemplate.content, true).firstElementChild;
      const wrapper = this.ownerDocument.importNode(this.#panelTemplate.content, true).firstElementChild;
      const id = `${this.#prefix}-${index}`, name = panel.getAttribute("name");
      panel.slot = id; wrapper.id = id; wrapper.querySelector("slot").name = id;
      button.id = `${id}-tab`; button.value = name; button.textContent = panel.getAttribute("label") || name;
      button.disabled = panel.hasAttribute("disabled"); button.setAttribute("aria-controls", id); wrapper.setAttribute("aria-labelledby", button.id);
      return { name, panel, button, wrapper };
    });
    this.#list.replaceChildren(this.#buttonTemplate, ...this.#tabs.map(({ button }) => button));
    this.#content.replaceChildren(this.#panelTemplate, ...this.#tabs.map(({ wrapper }) => wrapper));
    this.#sync();
  }

  #sync() {
    this.#list.setAttribute("aria-label", this.getAttribute("label") || "Sections");
    const active = this.#tabs.find((tab) => tab.name === this.value && !tab.button.disabled) || this.#tabs.find((tab) => !tab.button.disabled);
    if (active && active.name !== this.value) { this.value = active.name; return; }
    for (const tab of this.#tabs) {
      const selected = tab === active;
      tab.button.setAttribute("aria-selected", String(selected)); tab.button.tabIndex = selected ? 0 : -1;
      tab.panel.hidden = tab.wrapper.hidden = !selected;
    }
  }
}

function defineTabs() {
  if (!customElements.get("clinker-tabs")) customElements.define("clinker-tabs", Tabs);
}

const [markup, css] = await Promise.all([
  "<div data-editor-shell>\n<div class=\"cm-editor__resize-handle\" data-resize role=\"separator\" aria-label=\"Resize plot editor\" aria-orientation=\"horizontal\"></div>\n<header part=\"header\" class=\"cm-editor__header\">\n  <h2 class=\"cm-editor__title\">Plot editor</h2>\n  <div>\n    <button type=\"button\" data-command=\"undo\" disabled>Undo</button>\n    <button type=\"button\" data-command=\"redo\" disabled>Redo</button>\n  </div>\n</header>\n<clinker-tabs part=\"tabs\" data-editor-tabs label=\"Plot editor sections\" value=\"data\">\n  <clinker-data-panel part=\"data\" name=\"data\" label=\"Data\"></clinker-data-panel>\n  <clinker-appearance-panel part=\"appearance\" name=\"appearance\" label=\"Appearance\" hidden></clinker-appearance-panel>\n</clinker-tabs>\n</div>\n",
  ":host { display: flex; min-height: 180px; box-sizing: border-box; flex-direction: column; padding: 10px; border: 1px solid var(--cm-border,var(--_cm-border,#d7dce3)); border-radius: 10px; }\n:host([resizable]) { overflow: hidden; }\n.cm-editor__resize-handle { height: 10px; flex: 0 0 10px; margin: -10px -10px 0; cursor: ns-resize; touch-action: none; }\n.cm-editor__resize-handle::after { display: block; width: 36px; height: 3px; margin: 4px auto; border-radius: 3px; background: var(--cm-muted,var(--_cm-muted,#aaa)); content: \"\"; }\n.cm-editor__header { display: flex; align-items: center; justify-content: space-between; gap: 8px; margin-bottom: 8px; }\n.cm-editor__title { margin: 0; font-size: 17px; letter-spacing: -.01em; }\n[data-editor-shell] { display:flex; flex:1; flex-direction:column; min-width:0; min-height:0; }\n",
]);
const shellTemplate = createTemplate(markup);
const HTMLElementBase = globalThis.HTMLElement || class {};

/** Optional editor for an already rendered chart; edits use its public API. */
class ClinkerEditorElement extends HTMLElementBase {
  static observedAttributes = ["hidden"];
  #chart = null;
  #dataPanel = null;
  #appearancePanel = null;
  #unsubscribe = null;
  #resizeStart = null;

  constructor() {
    super();
    const root = this.attachShadow({ mode: "open" });
    const style = this.ownerDocument.createElement("style"); style.textContent = controlStyles + editorThemeStyles + css;
    root.append(style);
    root.addEventListener("click", (event) => {
      const command = event.target.closest("[data-command]")?.dataset.command;
      if (command === "undo") this.#chart?.undo?.();
      else if (command === "redo") this.#chart?.redo?.();
    });
    root.addEventListener("pointerdown", (event) => {
      if (!event.target.closest("[data-resize]") || !this.hasAttribute("resizable")) return;
      event.preventDefault();
      this.#resizeStart = { pointerY: event.clientY, height: this.getBoundingClientRect().height };
      event.target.setPointerCapture?.(event.pointerId);
    });
    this.addEventListener("pointermove", (event) => {
      if (!this.#resizeStart) return;
      const parentHeight = this.parentElement?.getBoundingClientRect().height || innerHeight;
      this.style.height = `${Math.max(180, Math.min(parentHeight - 48, this.#resizeStart.height - (event.clientY - this.#resizeStart.pointerY)))}px`;
    });
    this.addEventListener("pointerup", () => { this.#resizeStart = null; });
    this.addEventListener("pointercancel", () => { this.#resizeStart = null; });
    upgradeProperty(this, "chart");
  }

  connectedCallback() {
    this.classList.add("cm-ui", "cm-editor", "cm-data-ui");
    if (!this.#dataPanel) this.#renderShell();
    this.#subscribe(); this.#render();
  }

  disconnectedCallback() {
    this.#resizeStart = null;
    this.#unsubscribe?.(); this.#unsubscribe = null;
  }

  attributeChangedCallback(name, previous, next) {
    if (name !== "hidden" || previous === next || !this.isConnected) return;
    this.dispatchEvent(new CustomEvent(this.hidden ? "clinker-editor-close" : "clinker-editor-open", {
      bubbles: true, composed: true,
    }));
  }

  get chart() { return this.#chart; }
  set chart(value) {
    if (value === this.#chart) return;
    this.#unsubscribe?.(); this.#unsubscribe = null;
    this.#chart = value || null;
    this.#subscribe(); this.#render();
  }
  get value() { return this.#chart?.data?.(); }

  editAnnotations(panelUid) { this.#openData().editAnnotations(panelUid); return this; }
  editSelection(selection) { this.#openData().editSelection(selection); return this; }

  open({ tab = null } = {}) {
    if (tab !== null && !["data", "appearance"].includes(tab)) throw new TypeError(`Unknown editor tab: ${tab}`);
    if (!this.#dataPanel) { this.#renderShell(); this.#render(); }
    if (tab !== null) this.shadowRoot.querySelector("[data-editor-tabs]").value = tab;
    this.hidden = false;
    return this;
  }

  close() { this.hidden = true; return this; }

  async openAppearance({ section = null, focus = true } = {}) {
    if (!this.isConnected || !this.#chart) return false;
    if (!this.#appearancePanel.sections.includes(section) && section !== null) return false;
    this.open({ tab: "appearance" });
    return this.#appearancePanel.openSection(section, { focus });
  }

  #openData() {
    this.open({ tab: "data" });
    return this.#dataPanel;
  }

  reset() {
    this.#renderShell(); this.#render();
    return this;
  }

  #renderShell() {
    this.shadowRoot.querySelector("[data-editor-shell]")?.remove();
    this.shadowRoot.append(shellTemplate(this.ownerDocument));
    this.#dataPanel = this.shadowRoot.querySelector("clinker-data-panel");
    this.#appearancePanel = this.shadowRoot.querySelector("clinker-appearance-panel");
  }

  #subscribe() {
    if (this.#unsubscribe || !this.isConnected || !this.#chart?.on) return;
    this.#unsubscribe = this.#chart.on("change", (change) => {
      if (change.type === "theme.change") this.#syncTheme();
      if (change.type === "history.change") this.#renderHistory();
    });
  }

  #render() {
    if (!this.#dataPanel) return;
    this.#dataPanel.chart = this.#chart;
    this.#appearancePanel.chart = this.#chart;
    this.#syncTheme(); this.#renderHistory();
  }

  #renderHistory() {
    const undo = this.shadowRoot.querySelector('[data-command="undo"]');
    const redo = this.shadowRoot.querySelector('[data-command="redo"]');
    if (undo) undo.disabled = !this.#chart?.canUndo?.();
    if (redo) redo.disabled = !this.#chart?.canRedo?.();
  }

  #syncTheme() {
    applyEditorTheme(this, this.#chart);
  }

}

function defineClinkerEditor() {
  if (!globalThis.customElements) throw new Error("Custom elements are unavailable in this environment.");
  defineTabs(); defineDataPanel(); defineAppearancePanel();
  if (!customElements.get("clinker-editor")) customElements.define("clinker-editor", ClinkerEditorElement);
  return customElements.get("clinker-editor");
}

if (globalThis.customElements) defineClinkerEditor();

export { mountPlotSelectionToolbar };
