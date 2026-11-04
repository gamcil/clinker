import { i as installDataEditorStyle, a as installStyle, h as defineAnnotationDataForm, j as createAnnotationPanel, f as defineAnnotationEditor } from './editor/annotation-editor-3c63735d.mjs';
import { g as createActionRunner, e as colourHex, j as createControls, k as classifyControl } from './editor/search-field-595834d9.mjs';

/** Test world bounds after applying a preview translation. */

/** Convert a pointer event on a surface into chart-world coordinates. */
function worldPoint(surface, event, camera) {
  const bounds = surface.getBoundingClientRect();
  return {
    x: (event.clientX - bounds.left - camera.x) / camera.k,
    y: (event.clientY - bounds.top - camera.y) / camera.k,
  };
}

/** Quick styling stays in the plot; data and membership editing live in the data editor. */
function createAnnotationControls({ getData, patch, document, panelUids, annotationApi, onEditData }) {
  const { ANNOTATION_DEFAULTS: defaults, ANNOTATION_PALETTES: palettes, resolveAnnotationScale } = annotationApi;
  const node = document.createElement("fieldset"); node.className = "cm-annotation-controls";
  let applying = false;
  const run = createActionRunner((message) => { node.querySelector("[data-annotation-error]").textContent = message; });
  const current = () => (getData()?.annotationPanels || []).filter((panel) => panelUids.includes(panel.uid));
  const ui = createControls(document), { el, row } = ui;
  const input = (parent, label, field, value, type = "text", options = {}) => {
    const control = ui.input(row(parent, label), label, value, type, options);
    control.dataset.annotationField = field; return control;
  };
  const select = (parent, label, field, value, entries) => {
    const control = ui.select(row(parent, label), label, value, entries);
    control.dataset.annotationField = field; return control;
  };
  const paletteSelect = (parent, resolved) => {
    const name = Object.entries(palettes).find(([, colours]) => colours.length === resolved.colours.length &&
      colours.every((value, i) => colourHex(value) === colourHex(resolved.colours[i])))?.[0] || "custom";
    select(parent, "Palette", "palette", name, [["custom", "Custom"], ...Object.keys(palettes).map((key) => [key, key])]);
  };
  const gradient = (root, panel) => {
    const preview = root.querySelector("[data-colour-preview]");
    if (preview) preview.style.background = `linear-gradient(to right, ${resolveAnnotationScale(panel).colours.join(",")})`;
  };
  function build(panel, batch) {
    const node = document.createElement("fieldset");
    const groups = panel.type === "groups", bar = panel.type === "bar";
    input(node, "Title", "title", panel.title ?? panel.uid);
    const resolved = groups ? null : resolveAnnotationScale(panel);
    if (groups && batch) {
      const colours = panel.groups.map((group) => colourHex(group.colour ?? panel.colour ?? defaults.colour));
      const control = input(node, "Colour", "colour", colours[0], "color");
      if (colours.some((colour) => colour !== colours[0])) control.dataset.mixed = "true";
    }
    else if (groups) panel.groups.forEach((group, i) => input(node, group.label || group.uid, `group:${i}:colour`,
      group.colour ?? panel.colour ?? defaults.colour, "color"));
    else if (bar) input(node, "Colour", "colour", panel.colour ?? defaults.colour, "color");
    else if (resolved.type === "categorical") resolved.domain.forEach((name, i) => input(node, name, `colour:${i}`, resolved.colour(name), "color"));
    else {
      paletteSelect(node, resolved);
      const preview = el("div", undefined, node); preview.dataset.colourPreview = ""; preview.className = "cm-annotation-gradient";
      preview.setAttribute("role", "img"); preview.setAttribute("aria-label", "Colour gradient"); gradient(node, panel);
    }
    if (resolved?.type === "continuous") {
      const range = el("div", undefined, node); range.className = "cm-annotation-range";
      select(range, "Range", "autoDomain", panel.scale?.domain ? "custom" : "auto", [["auto", "Auto"], ["custom", "Custom"]]);
      if (panel.scale?.domain) {
        input(range, "Minimum", "domainMin", resolved.domain[0], "number");
        input(range, "Maximum", "domainMax", resolved.domain[1], "number");
      }
    }
    if (!groups) input(node, bar ? "Show axis" : "Show legend", "showKey", panel.showKey !== false, "checkbox");
    const data = el("button", "Edit data…", node); data.type = "button"; data.onclick = () => onEditData?.();
    const advanced = el("details", undefined, node); el("summary", "Advanced", advanced);
    const sizes = groups ? [["markerWidth", "Bar thickness", 3], ["laneGap", "Lane gap", 12]]
      : bar ? [["width", "Panel width", defaults.width], ["cellHeight", "Bar height", defaults.cellHeight]]
        : [["cellWidth", "Cell width", defaults.cellWidth], ["cellHeight", "Cell height", defaults.cellHeight], ["cellGap", "Cell gap", defaults.cellGap]];
    const fonts = [["titleFontSize", "Title font size", defaults.titleFontSize],
      ...(!groups ? [["columnFontSize", "Column title font size", defaults.columnFontSize]] : []),
      ["labelFontSize", groups ? "Group label font size" : "Axis/legend font size", defaults.labelFontSize]];
    for (const [field, label, fallback] of [...sizes, ["gap", "Panel gap", defaults.gap], ...fonts]) {
      input(advanced, label, field, panel[field] ?? fallback, "number", { min: ["gap", "cellGap"].includes(field) ? "0" : "0.1" });
    }
    if (bar) input(advanced, "Show y-axis line", "showYAxisLine", panel.showYAxisLine ?? false, "checkbox");
    if (bar) {
      const spacing = input(advanced, "X-axis tick spacing", "tickStep", panel.tickStep ?? "", "number", { min: "0" });
      spacing.placeholder = "Auto"; spacing.title = "Leave empty for automatic tick spacing.";
    }
    if (bar) input(advanced, "Show x-axis grid lines", "showXGridLines", panel.showXGridLines ?? false, "checkbox");
    if (!groups) {
      input(advanced, "Border width", "cellBorderWidth", panel.cellBorderWidth ?? defaults.cellBorderWidth, "number", { min: "0" });
      input(advanced, "Border colour", "cellBorderColour", panel.cellBorderColour ?? defaults.cellBorderColour, "color");
      input(advanced, "Missing colour", "missingColour", panel.missingColour ?? defaults.missingColour, "color");
      if (!bar && resolved.type === "continuous") resolved.colours.forEach((value, i) => input(advanced, `Colour stop ${i + 1}`, `colour:${i}`, value, "color"));
      else if (!bar) paletteSelect(advanced, resolved);
      input(advanced, "Show column labels", "showColumnLabels", panel.showColumnLabels !== false, "checkbox");
      input(advanced, "Heading angle", "headingRotation", panel.headingRotation ?? (bar ? 0 : defaults.headingRotation), "number", { min: "-90", max: "0" });
    }
    return node;
  }
  function render() {
    const open = node.querySelector("details")?.open || false, panels = current();
    node.replaceChildren(); if (!panels.length) return;
    const forms = panels.map((panel) => build(panel, panels.length > 1)), first = forms[0];
    if (panels.length > 1) {
      first.querySelector("button")?.remove();
      first.querySelector("[data-colour-preview]")?.remove();
      const fields = forms.map((form) => new Map([...form.querySelectorAll("[data-annotation-field]")]
        .map((control) => [control.dataset.annotationField, control])));
      for (const [field, control] of fields[0]) {
        const peers = fields.map((map) => map.get(field));
        const label = { cellHeight: "Row height", labelFontSize: "Label font size", showKey: "Show axis/legend" }[field];
        if (field === "title" || field.startsWith("group:") || peers.some((peer) => !peer ||
            !label && peer.getAttribute("aria-label") !== control.getAttribute("aria-label"))) {
          control.closest("label").remove(); continue;
        }
        if (label && peers.some((peer) => peer.getAttribute("aria-label") !== control.getAttribute("aria-label"))) {
          control.setAttribute("aria-label", label); control.closest("label").querySelector("span").textContent = label;
        }
        const value = (input) => input.type === "checkbox" ? input.checked : input.value;
        if (peers.every((peer) => !peer.dataset.mixed && value(peer) === value(control))) continue;
        control.dataset.mixed = "true";
        if (control.type === "checkbox") { control.checked = false; control.indeterminate = true; control.setAttribute("aria-checked", "mixed"); }
        else if (control.tagName === "SELECT") {
          const option = el("option", "Mixed", control); option.value = ""; option.disabled = true; control.value = "";
        } else if (control.type === "color") el("small", "Mixed", control.closest("label")).dataset.mixedLabel = "";
        else { control.value = ""; control.placeholder = "Mixed"; }
      }
      first.querySelectorAll(".cm-annotation-range, details").forEach((section) => {
        if (!section.querySelector("[data-annotation-field]")) section.remove();
      });
    }
    node.append(...first.childNodes);
    if (node.querySelector("details")) node.querySelector("details").open = open;
    const error = el("div", undefined, node); error.dataset.annotationError = ""; error.setAttribute("role", "alert");
  }
  const edit = (event) => {
    const control = event.target, field = control.dataset.annotationField;
    if (!field || (control.type === "text") !== (event.type === "change")) return;
    control.setCustomValidity("");
    if (!control.checkValidity() || control.type === "number" && control.value === "" && field !== "tickStep") return;
    const panels = current();
    if (!panels.length) return;
    const value = control.type === "checkbox" ? control.checked : control.type === "number" ? (control.value === "" ? null : control.valueAsNumber) : control.value;
    let refresh = false;
    const operations = panels.map((panel) => {
      const resolved = panel.type === "groups" ? null : resolveAnnotationScale(panel);
      let changes;
      const scale = { ...panel.scale };
      if (field.startsWith("group:")) {
        const index = Number(field.split(":")[1]);
        changes = { groups: panel.groups.map((group, i) => i === index ? { ...group, colour: value } : group) };
      } else if (field === "autoDomain") {
        if (value === "auto") delete scale.domain; else scale.domain = resolved.domain;
        changes = { scale }; refresh = true;
      } else if (field === "domainMin" || field === "domainMax") {
        const domain = [...resolved.domain]; domain[field === "domainMin" ? 0 : 1] = value;
        changes = { scale: { ...scale, domain } };
      } else if (field === "palette") {
        if (value === "custom") return null;
        changes = { scale: { ...scale, colours: [...palettes[value]] } }; refresh = true;
      } else if (field.startsWith("colour:")) {
        const colours = resolved.type === "categorical" ? resolved.domain.map(resolved.colour) : [...resolved.colours];
        colours[Number(field.split(":")[1])] = value; changes = { scale: { ...scale, colours } };
      } else if (field === "colour" && panel.type === "groups") {
        changes = { colour: value, groups: panel.groups.map((group) => ({ ...group, colour: value })) };
      } else changes = { [field]: value };
      return { type: "annotationPanels.update", ids: [panel.uid], changes };
    }).filter(Boolean);
    if (!operations.length) return;
    applying = true;
    try { run(() => {
      patch(operations);
      delete control.dataset.mixed; control.removeAttribute("aria-checked");
      control.closest("label")?.querySelector("[data-mixed-label]")?.remove();
      node.querySelectorAll("input").forEach((input) => input.setCustomValidity(""));
      if (refresh) render();
      else { gradient(node, current()[0]); if (field.startsWith("colour:")) { const palette = node.querySelector('[data-annotation-field="palette"]'); if (palette) palette.value = "custom"; } }
    }, control); } finally { applying = false; }
  };
  node.addEventListener("input", edit); node.addEventListener("change", edit);
  render();
  return { node, update() { if (!applying) render(); } };
}

const css = ".cm-annotation-ui { position:fixed; z-index:4; box-sizing:border-box; --_cm-control-height:24px; --_cm-control-padding:2px 5px; --_cm-control-radius:4px; color:var(--cm-foreground,var(--_cm-foreground,#273444)); font-size:12px; line-height:1.3; }\n.cm-annotation-selected { z-index:3; }\n.cm-annotation-ui[hidden] { display:none!important; }\n.cm-annotation-frame, .cm-annotation-selected { pointer-events:none; border:1px solid var(--cm-accent,var(--_cm-accent)); border-radius:4px; background:rgba(22,119,255,.025); }\n.cm-annotation-frame button { position:absolute; right:0; top:-28px; pointer-events:auto; }\n.cm-annotation-frame:not(.cm-annotation-add)::before { content:\"\"; position:absolute; top:-28px; right:0; width:48px; height:28px; pointer-events:auto; }\n.cm-annotation-add { border-style:dashed; background:transparent; }\n.cm-annotation-add button { top:50%; right:50%; transform:translate(50%,-50%); font-size:22px; width:32px; height:32px; }\n.cm-annotation-actions { display:flex; gap:5px; border-top:1px solid var(--cm-border,var(--_cm-border,#e5e9ef)); margin-top:8px; padding-top:8px; }\n.cm-annotation-menu p { font-size:11px; color:var(--cm-muted,var(--_cm-muted,#607086)); }\n.cm-annotation-menu { z-index:5; width:280px; max-width:calc(100vw - 16px); max-height:calc(100vh - 16px); overflow:auto; padding:8px; border:1px solid var(--cm-border,var(--_cm-border,#cbd3de)); border-radius:4px; background:var(--cm-surface,var(--_cm-surface,#fff)); box-shadow:0 2px 8px rgba(0,0,0,.2); }\n.cm-annotation-menu header { display:flex; align-items:center; gap:8px; margin-bottom:6px; }\n.cm-annotation-menu header strong { flex:1; min-width:0; overflow-wrap:anywhere; }\n.cm-annotation-menu header button { min-width:24px; padding:0; }\n.cm-annotation-menu fieldset { margin:0; padding:0; border:0; min-width:0; }\n.cm-annotation-menu legend { display:none; }\n.cm-annotation-menu label { display:flex; gap:6px; align-items:center; margin:5px 0; }\n.cm-annotation-menu label > span { flex:1; min-width:0; overflow-wrap:anywhere; }\n.cm-annotation-menu input[type=number] { width:72px; }\n.cm-annotation-menu input[type=text],.cm-annotation-menu select { width:132px; }\n.cm-annotation-menu input[type=color] { width:28px; height:22px; padding:0; border:0; background:transparent; }\n.cm-annotation-menu input[type=checkbox] { width:13px; height:13px; margin:0; }\n.cm-annotation-menu details:not(.cm-annotation-data-editor *) { border-top:1px solid var(--cm-border,var(--_cm-border,#e5e9ef)); margin-top:8px; }\n.cm-annotation-menu summary:not(.cm-annotation-data-editor *) { color:var(--cm-muted,var(--_cm-muted,#4f6d95)); font-weight:600; padding-top:6px; cursor:pointer; }\n.cm-annotation-menu [data-annotation-error] { color:var(--cm-danger,var(--_cm-danger,#b42318)); }\n.cm-annotation-gradient { height:12px; border:1px solid var(--cm-border,var(--_cm-border,#cbd3de)); border-radius:3px; margin:6px 0; }\n.cm-annotation-range { display:flex; gap:6px; align-items:center; }\n.cm-annotation-range label:first-child { flex:1; }\n.cm-annotation-range label:first-child select { width:78px; }\n.cm-annotation-range label:not(:first-child) > span { display:none; }\n.cm-annotation-range input[type=number] { width:60px; }\n.cm-annotation-menu--data { --_cm-control-height:30px; --_cm-control-padding:4px 8px; --_cm-control-radius:6px; width:760px; display:flex; flex-direction:column; height:min(600px,calc(100vh - 16px)); box-sizing:border-box; }\n.cm-annotation-menu--data > .cm-annotation-data-editor { flex:1; min-height:0; }\n.cm-annotation-hidden { border:1px solid var(--cm-border,var(--_cm-border,#cbd3de)); border-radius:4px; padding:3px 6px; background:var(--cm-surface,var(--_cm-surface,#fff)); cursor:pointer; }\n.cm-annotation-menu .cm-annotation-restore { display:block; width:100%; margin:5px 0; text-align:left; }\n.annotation-tooltip { position:fixed; pointer-events:none; z-index:6; max-width:300px; padding:8px 10px; border:1px solid var(--cm-border,var(--_cm-border)); border-radius:4px; background:var(--cm-surface,var(--_cm-surface)); box-shadow:0 2px 8px rgba(0,0,0,.15); color:var(--cm-foreground,var(--_cm-foreground)); font-size:12px; line-height:1.5; overflow-wrap:anywhere; }\n\n.cm-annotation-selection { margin-bottom:8px; }\n.cm-annotation-selection > label { margin-left:6px; }\n.cm-annotation-controls [data-mixed-label] { color:var(--cm-muted,var(--_cm-muted)); }\n";

/** Renderer-independent editing affordances, kept outside the exported scene. */
function createAnnotationOverlay({ container, getData, patch, annotationApi, editData }) {
  const getPanels = () => getData()?.annotationPanels || [];
  const { annotationPanelsForFrame, annotationPanelBounds, hitTestAnnotationCell } = annotationApi;
  const document = container.ownerDocument, window = document.defaultView;
  const styleRoot = container.getRootNode();
  installDataEditorStyle(styleRoot);
  installStyle(styleRoot, "clustermap-annotation-ui-style", css);
  const element = (tag, className, parent = container) => {
    const node = document.createElement(tag); node.className = className; classifyControl(node); parent.append(node); return node;
  };
  const frame = element("div", "cm-ui cm-data-ui cm-annotation-ui cm-annotation-frame"); frame.hidden = true;
  const edit = element("button", "", frame); edit.type = "button"; edit.textContent = "Edit";
  edit.setAttribute("aria-haspopup", "dialog"); edit.setAttribute("aria-expanded", "false");
  const addFrame = element("div", "cm-ui cm-data-ui cm-annotation-ui cm-annotation-frame cm-annotation-add"); addFrame.hidden = true;
  const add = element("button", "", addFrame); add.type = "button"; add.textContent = "+";
  add.setAttribute("aria-label", "Add annotation panel"); add.setAttribute("aria-haspopup", "dialog");
  const hiddenButton = element("button", "cm-ui cm-data-ui cm-annotation-ui cm-annotation-hidden");
  hiddenButton.type = "button"; hiddenButton.hidden = true; hiddenButton.setAttribute("aria-haspopup", "dialog");
  const menu = element("div", "cm-ui cm-data-ui cm-annotation-ui cm-annotation-menu"); menu.hidden = true;
  menu.setAttribute("role", "dialog"); menu.setAttribute("aria-label", "Annotation panel settings");
  const tooltip = element("div", "cm-ui annotation-tooltip"); tooltip.hidden = true;
  tooltip.setAttribute("role", "tooltip");
  let context = null, surfaces = [], hovered = null, menuState = null;
  let tooltipKey = null, selected = new Set();
  const selectionFrames = new Map();
  let panelAreas = [], addArea = null;
  let areaScene = null, areaPreview = null;
  const contains = (bounds, point) => point.x >= bounds.minX && point.x <= bounds.maxX && point.y >= bounds.minY && point.y <= bounds.maxY;
  const hideFrame = () => { if (!menuState) { hovered = null; frame.hidden = true; addFrame.hidden = true; } };
  const hideTooltip = () => { tooltip.hidden = true; tooltipKey = null; };
  const updateTooltip = (event, point) => {
    const cell = hitTestAnnotationCell(context.getScene(), point, context.getPreview());
    if (!cell) { hideTooltip(); return; }
    const key = JSON.stringify([cell.panelUid, cell.clusterUid, cell.column, cell.value]);
    if (tooltipKey !== key) {
      tooltipKey = key; tooltip.replaceChildren();
      element("strong", "", tooltip).textContent = cell.title || cell.panelUid;
      for (const value of [cell.clusterName, `${cell.column}: ${cell.value === null ? "Missing" : cell.value}`]) {
        element("div", "", tooltip).textContent = value;
      }
    }
    tooltip.hidden = false;
    const bounds = tooltip.getBoundingClientRect();
    tooltip.style.left = `${Math.max(8, Math.min(event.clientX + 12, window.innerWidth - bounds.width - 8))}px`;
    tooltip.style.top = `${Math.max(8, Math.min(event.clientY + 12, window.innerHeight - bounds.height - 8))}px`;
  };
  const positionMenu = () => {
    if (!menuState?.anchor) return;
    const { anchor } = menuState;
    const bounds = menu.getBoundingClientRect();
    const right = anchor.right + 8;
    const x = right + bounds.width <= window.innerWidth - 8 ? right : anchor.left - bounds.width - 8;
    menu.style.left = `${Math.max(8, Math.min(x, window.innerWidth - bounds.width - 8))}px`;
    menu.style.top = `${Math.max(8, Math.min(anchor.top, window.innerHeight - bounds.height - 8))}px`;
  };
  const close = (returnFocus = false) => {
    hideTooltip();
    if (!menuState && !hovered && frame.hidden && addFrame.hidden) return;
    const opener = menuState?.opener;
    menuState = null; hovered = null; menu.hidden = true; frame.hidden = true; addFrame.hidden = true;
    menu.replaceChildren();
    for (const button of [edit, hiddenButton, add]) button.setAttribute("aria-expanded", "false");
    if (returnFocus && opener?.isConnected && !opener.hidden && !opener.closest("[hidden]")) opener.focus({ preventScroll: true });
    else if (returnFocus) surfaces[0]?.focus?.({ preventScroll: true });
  };
  const syncDataEditor = (state) => {
    const data = getData();
    Object.assign(state, { panels: data.annotationPanels, panel: data.annotationPanels?.find((panel) => panel.uid === state.editor.value), clusters: data.clusters });
    state.editor.state = { data, clusters: data.clusters, annotationApi,
      patch: (operations) => { patch(operations); syncDataEditor(state); } };
  };
  const beginMenu = (kind, title, button, uid = null) => {
    hideTooltip();
    // Keep the menu anchored when creating a panel moves its opener offscreen.
    const anchor = button.closest("[hidden]") ? menuState?.anchor : button.getBoundingClientRect();
    menuState = { kind, uid, anchor, opener: button };
    hovered = null; frame.hidden = true; addFrame.hidden = true;
    for (const opener of [edit, hiddenButton, add]) opener.setAttribute("aria-expanded", String(opener === button));
    menu.replaceChildren();
    menu.classList.remove("cm-annotation-menu--data");
    const row = element("header", "", menu), heading = element("strong", "", row);
    heading.textContent = title;
    const closeButton = element("button", "", row); closeButton.type = "button";
    closeButton.textContent = "×"; closeButton.setAttribute("aria-label", "Close annotation settings");
    closeButton.onclick = () => close(true);
    menu.hidden = false;
    return menuState;
  };
  const screenArea = (area, surface) => {
    const camera = context.getCamera(), rect = surface.getBoundingClientRect();
    const left = Math.max(rect.left + 2, rect.left + camera.x + area.minX * camera.k - 10, 0);
    const top = Math.max(rect.top + 30, rect.top + camera.y + area.minY * camera.k - 12, 30);
    const right = Math.min(rect.right - 2, rect.left + camera.x + area.maxX * camera.k + 10, window.innerWidth);
    const bottom = Math.min(rect.bottom - 2, rect.top + camera.y + area.maxY * camera.k + 12, window.innerHeight);
    return right > left && bottom > top ? { left, top, right, bottom, width: right - left, height: bottom - top } : null;
  };
  const placeFrame = (node, rect) => {
    node.hidden = !rect;
    if (rect) Object.assign(node.style, { left: `${rect.left}px`, top: `${rect.top}px`, width: `${rect.width}px`, height: `${rect.height}px` });
  };
  const showAddFrame = () => placeFrame(addFrame, addScreenArea());
  const addScreenArea = () => {
    if (!addArea || !surfaces[0]) return null;
    const surface = surfaces[0].getBoundingClientRect(), camera = context.getCamera();
    const edge = surface.left + camera.x + (addArea.minX - 24) * camera.k;
    // A fitted figure can leave less room than the empty slot needs. Keep a
    // small hover target at its visible right edge without changing export bounds.
    const right = Math.min(surface.right - 2, window.innerWidth - 2);
    if (edge > right || edge < surface.left) return screenArea(addArea, surfaces[0]);
    const rect = screenArea(addArea, surfaces[0]);
    if (rect && rect.width >= 44) return rect;
    const vertical = screenArea({ ...addArea, minX: (right - 44 - surface.left - camera.x) / camera.k,
      maxX: (right - surface.left - camera.x) / camera.k }, surfaces[0]);
    return vertical ? { ...vertical, left: Math.max(surface.left + 2, right - 44), right,
      width: Math.min(44, right - surface.left - 2) } : null;
  };
  const showFrame = (uid) => {
    const area = panelAreas.find((area) => area.uid === uid);
    const rect = area && surfaces[0] && screenArea(area.bounds, surfaces[0]);
    placeFrame(frame, rect);
    if (!rect) return;
    edit.setAttribute("aria-label", selected.has(uid) && selected.size > 1
      ? `Edit ${selected.size} annotation panels` : `Edit annotation panel ${area.title || uid}`);
  };
  const openPanel = (uid, button) => {
    const panel = getPanels().find((panel) => panel.uid === uid);
    if (!panel) return;
    if (!selected.has(uid)) selected = new Set([uid]);
    const ids = getPanels().filter((panel) => selected.has(panel.uid)).map((panel) => panel.uid);
    const state = beginMenu("panel", ids.length > 1 ? `${ids.length} annotation panels` : panel.title || uid, button, uid);
    state.ids = ids; state.panelOrder = getPanels().map((panel) => panel.uid);
    if (getPanels().length > 1) {
      const chooser = element("details", "cm-annotation-selection", menu);
      element("summary", "", chooser).textContent = `Panels (${ids.length} selected)`;
      const selectAll = element("button", "", chooser); selectAll.type = "button"; selectAll.textContent = "Select all panels";
      const choose = (next, focusUid) => {
        selected = new Set(next); openPanel(next[0], button);
        menu.querySelector(".cm-annotation-selection").open = true;
        [...menu.querySelectorAll("[data-panel-select]")].find((input) => input.dataset.panelSelect === focusUid)?.focus({ preventScroll: true });
      };
      selectAll.onclick = () => choose(getPanels().map((panel) => panel.uid), uid);
      for (const item of getPanels()) {
        const label = element("label", "", chooser);
        element("span", "", label).textContent = `${item.title || item.uid}${item.show === false ? " (hidden)" : ""}`;
        const check = element("input", "", label); check.type = "checkbox";
        check.dataset.panelSelect = item.uid; check.checked = selected.has(item.uid);
        check.disabled = check.checked && ids.length === 1;
        check.setAttribute("aria-label", `Select annotation panel ${item.title || item.uid}`);
        check.onchange = () => choose(check.checked ? [...ids, item.uid] : ids.filter((id) => id !== item.uid), item.uid);
      }
    }
    state.controls = createAnnotationControls({ getData, patch, document, panelUids: ids, annotationApi,
      onEditData: () => {
        if (editData?.(uid)) { close(); return; }
        const dataState = beginMenu("data", "Annotation data", edit, uid);
        menu.classList.add("cm-annotation-menu--data");
        defineAnnotationEditor();
        dataState.editor = document.createElement("clinker-annotation-editor");
        syncDataEditor(dataState);
        dataState.editor.select(uid); menu.append(dataState.editor); positionMenu();
      } });
    menu.append(state.controls.node);
    const actions = element("div", "cm-annotation-actions", menu);
    const movePanel = (direction) => {
      const ids = getPanels().map((panel) => panel.uid), index = ids.indexOf(uid), target = index + direction;
      if (index < 0 || target < 0 || target >= ids.length) return;
      [ids[index], ids[target]] = [ids[target], ids[index]];
      patch([{ type: "annotationPanels.reorder", ids }]);
    };
    const action = (label, callback) => {
      const button = element("button", "", actions); button.type = "button"; button.textContent = label; button.onclick = callback; return button;
    };
    state.actions = { left: action("←", () => movePanel(-1)), right: action("→", () => movePanel(1)),
      visibility: action(panel.show === false ? "Show" : "Hide", () => patch([{ type: "annotationPanels.update", ids,
        changes: { show: ids.every((id) => getPanels().find((item) => item.uid === id).show === false) } }])) };
    state.actions.left.setAttribute("aria-label", "Move panel left"); state.actions.left.title = "Move panel left";
    state.actions.right.setAttribute("aria-label", "Move panel right"); state.actions.right.title = "Move panel right";
    action("Remove", () => patch([{ type: "annotationPanels.delete", ids }]));
    refresh();
    positionMenu();
    state.controls.node.querySelector("input, select")?.focus({ preventScroll: true });
  };
  edit.onclick = () => { if (hovered) openPanel(hovered, edit); };
  add.onclick = () => {
    beginMenu("add", "Add annotation panel", add);
    defineAnnotationDataForm();
    const form = document.createElement("clinker-annotation-data-form"); form.clusters = getData().clusters;
    form.addEventListener("data-submit", (event) => {
      try {
        const panel = createAnnotationPanel(event.detail, getPanels());
        patch([{ type: "annotationPanels.create", panel }]); openPanel(panel.uid, add);
      } catch (exception) { form.error = exception.message; }
    });
    menu.append(form); showAddFrame(); positionMenu();
    form.focus({ preventScroll: true });
  };
  hiddenButton.onclick = () => {
    beginMenu("hidden", "Hidden annotation panels", hiddenButton);
    for (const panel of getPanels().filter((panel) => panel.show === false)) {
      const button = element("button", "cm-annotation-restore", menu); button.type = "button";
      button.textContent = `Edit ${panel.title || panel.uid}`;
      button.onclick = () => openPanel(panel.uid, hiddenButton);
    }
    positionMenu();
    menu.querySelector("button")?.focus({ preventScroll: true });
  };
  const refresh = (event) => {
    if (event) hideTooltip();
    if (!context) return;
    const scene = context.getScene();
    const preview = context.getPreview();
    if (scene !== areaScene || preview !== areaPreview) {
      areaScene = scene; areaPreview = preview;
      panelAreas = scene ? annotationPanelsForFrame(scene, preview).map((panel) =>
        ({ uid: panel.uid, title: panel.title, bounds: annotationPanelBounds(panel) })) : [];
      const bounds = scene?.figureBounds || scene?.bounds;
      addArea = bounds ? { minX: bounds.maxX + 24, maxX: bounds.maxX + 120,
        minY: panelAreas.length ? Math.min(...panelAreas.map((area) => area.bounds.minY)) : bounds.minY,
        maxY: panelAreas.length ? Math.max(...panelAreas.map((area) => area.bounds.maxY)) : bounds.maxY } : null;
    }
    const panels = getPanels();
    if (menuState?.kind === "panel") {
      const choices = [...menu.querySelectorAll("[data-panel-select]")];
      for (const check of choices) {
        const panel = panels.find((panel) => panel.uid === check.dataset.panelSelect);
        if (panel) check.previousElementSibling.textContent = `${panel.title || panel.uid}${panel.show === false ? " (hidden)" : ""}`;
      }
    }
    selected = new Set([...selected].filter((uid) => panels.some((panel) => panel.uid === uid)));
    for (const [uid, outline] of selectionFrames) if (!selected.has(uid)) { outline.remove(); selectionFrames.delete(uid); }
    for (const uid of selected) {
      let outline = selectionFrames.get(uid);
      if (!outline) { outline = element("div", "cm-ui cm-annotation-ui cm-annotation-selected"); selectionFrames.set(uid, outline); }
      const area = panelAreas.find((area) => area.uid === uid);
      const rect = area && surfaces[0] && screenArea(area.bounds, surfaces[0]);
      placeFrame(outline, rect);
    }
    if (menuState?.kind === "panel" && (menuState.ids.some((uid) => !selected.has(uid)) ||
        menuState.panelOrder.length !== panels.length || panels.some((panel, i) => panel.uid !== menuState.panelOrder[i]))) {
      if (selected.size) {
        const choosing = menu.querySelector(".cm-annotation-selection")?.open;
        openPanel([...selected][0], menuState.opener);
        if (choosing && menu.querySelector(".cm-annotation-selection")) menu.querySelector(".cm-annotation-selection").open = true;
      } else close();
    }
    const active = menuState?.kind === "panel" ? menuState.uid : null;
    if (active || hovered) showFrame(active || hovered);
    if (menuState?.kind === "add") showAddFrame();
    if (active) {
      const ids = menuState.ids, batch = ids.length > 1;
      menu.querySelector("header strong").textContent = batch ? `${ids.length} annotation panels` : panels.find((panel) => panel.uid === active)?.title || active;
      menuState.actions.left.hidden = menuState.actions.right.hidden = batch;
      const index = panels.findIndex((panel) => panel.uid === active);
      menuState.actions.left.disabled = index <= 0; menuState.actions.right.disabled = index >= panels.length - 1;
      menuState.actions.visibility.textContent = ids.every((uid) => panels.find((panel) => panel.uid === uid).show === false) ? "Show" : "Hide";
    }
    const hidden = panels.filter((panel) => panel.show === false);
    hiddenButton.hidden = !hidden.length || !surfaces[0];
    hiddenButton.textContent = `Hidden panels (${hidden.length})`;
    if (!hiddenButton.hidden && scene?.bounds) {
      const surface = surfaces[0].getBoundingClientRect(), camera = context.getCamera();
      const top = panelAreas.length ? Math.min(...panelAreas.map((area) => area.bounds.minY)) : scene.bounds.minY;
      const x = surface.left + camera.x + ((scene.annotationAnchorX ?? scene.bounds.maxX) + 24) * camera.k;
      const y = surface.top + camera.y + top * camera.k - 30;
      hiddenButton.style.left = `${Math.max(surface.left + 8, Math.min(x, surface.right - hiddenButton.offsetWidth - 8))}px`;
      hiddenButton.style.top = `${Math.max(surface.top + 8, Math.min(y, surface.bottom - 32))}px`;
    }
    positionMenu();
  };
  const move = (event) => {
    if (menuState) { hideTooltip(); return; }
    if (event.buttons || context.isBusy()) { hideTooltip(); hideFrame(); return; }
    refresh();
    const camera = context.getCamera();
    const point = worldPoint(event.currentTarget, event, camera);
    updateTooltip(event, point);
    const padded = (bounds) => ({ minX: bounds.minX - 10 / camera.k, maxX: bounds.maxX + 10 / camera.k,
      minY: bounds.minY - 12 / camera.k, maxY: bounds.maxY + 12 / camera.k });
    const area = [...panelAreas].reverse().find((area) => contains(padded(area.bounds), point));
    hovered = area?.uid || null;
    if (hovered) { showFrame(hovered); addFrame.hidden = true; }
    else {
      const slot = addScreenArea();
      if (slot && event.clientX >= slot.left && event.clientX <= slot.right && event.clientY >= slot.top && event.clientY <= slot.bottom) {
        frame.hidden = true; showAddFrame();
      } else hideFrame();
    }
  };
  const leave = (event) => {
    hideTooltip();
    if (frame.contains(event.relatedTarget) || addFrame.contains(event.relatedTarget) || menu.contains(event.relatedTarget)) return;
    hideFrame();
  };
  const outside = (event) => {
    const path = event.composedPath();
    if (menuState && !path.includes(menu) && !path.includes(frame) && !path.includes(addFrame) && event.target !== hiddenButton) close();
  };
  const dismissForGesture = (event) => {
    if (event.type === "pointerdown" && !event.button && (event.ctrlKey || event.metaKey || event.shiftKey) && !context.isBusy()) {
      refresh();
      const point = worldPoint(event.currentTarget, event, context.getCamera());
      const area = [...panelAreas].reverse().find((area) => contains(area.bounds, point));
      if (area) {
        if (selected.has(area.uid)) selected.delete(area.uid); else selected.add(area.uid);
        event.preventDefault(); event.stopImmediatePropagation(); close(); hovered = area.uid; refresh();
        return;
      }
    }
    if (event.type === "pointerdown") selected.clear();
    if (menu.contains(document.activeElement)) document.activeElement.blur();
    close(); refresh();
  };
  const escape = (event) => {
    if (event.key !== "Escape" || event.defaultPrevented) return;
    if (!menuState && !selected.size && !event.composedPath().includes(container) && !container.contains(document.activeElement)) return;
    const open = Boolean(menuState);
    if (!open) selected.clear();
    close(open); refresh();
    if (open) event.preventDefault();
  };
  const detach = () => {
    for (const surface of surfaces) {
      surface.removeEventListener("pointermove", move); surface.removeEventListener("pointerleave", leave);
      surface.removeEventListener("pointerdown", dismissForGesture, true); surface.removeEventListener("wheel", dismissForGesture);
    }
  };
  frame.addEventListener("pointerleave", leave);
  addFrame.addEventListener("pointerleave", leave);
  menu.addEventListener("toggle", positionMenu, true);
  window.addEventListener("click", outside); window.addEventListener("keydown", escape);
  window.addEventListener("resize", refresh); window.addEventListener("scroll", refresh, true);
  return {
    get isOpen() { return Boolean(menuState); },
    refresh,
    update(next) {
      const nextSurfaces = next.surfaces.filter(Boolean);
      const surfacesChanged = surfaces.length !== nextSurfaces.length || surfaces.some((surface, i) => surface !== nextSurfaces[i]);
      if (surfacesChanged) detach();
      context = next; surfaces = nextSurfaces;
      if (surfacesChanged) for (const surface of surfaces) {
        surface.addEventListener("pointermove", move); surface.addEventListener("pointerleave", leave);
        surface.addEventListener("pointerdown", dismissForGesture, true); surface.addEventListener("wheel", dismissForGesture);
      }
      hideTooltip();
      if (menuState?.editor) {
        const data = getData(), panel = data.annotationPanels?.find((panel) => panel.uid === menuState.editor.value);
        if (menuState.panels !== data.annotationPanels || menuState.panel !== panel || menuState.clusters !== data.clusters) {
          syncDataEditor(menuState);
        }
      } else menuState?.controls?.update();
      refresh();
    },
    dispose() {
      detach(); window.removeEventListener("click", outside); window.removeEventListener("keydown", escape);
      window.removeEventListener("resize", refresh); window.removeEventListener("scroll", refresh, true);
      frame.remove(); addFrame.remove(); menu.remove(); hiddenButton.remove(); tooltip.remove();
      for (const outline of selectionFrames.values()) outline.remove();
    },
  };
}

export { createAnnotationOverlay };
