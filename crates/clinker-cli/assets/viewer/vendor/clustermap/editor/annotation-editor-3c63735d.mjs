import { b as uiAppearanceStyles, a as controlStyles, l as lightControlStyles, m as mountComponent, c as createTemplate, d as updateSelectOptions, e as colourHex, f as defineSearchField, g as createActionRunner } from './search-field-595834d9.mjs';

/** Rollup embeds these files; unbundled browser modules load them by URL. */

function installStyle(root, id, css) {
  if (root.querySelector(`style[id="${id}"]`)) return;
  const style = (root.ownerDocument || root).createElement("style");
  style.id = id; style.textContent = css;
  (root.head || root).append(style);
}

/** Shared fixed-height row window for scrollable data editors. */
const DATA_ROW_HEIGHT = 34;

function virtualRowWindow(count, scrollTop, height, { rowHeight = DATA_ROW_HEIGHT, overscan = 8, threshold = 200 } = {}) {
  const virtual = count > threshold;
  const first = virtual ? Math.min(Math.max(0, count - 1), Math.max(0, Math.floor(scrollTop / rowHeight) - overscan)) : 0;
  const last = virtual ? Math.min(count, first + Math.ceil((height || 320) / rowHeight) + overscan * 2) : count;
  return { first, last, before: first * rowHeight, after: (count - last) * rowHeight };
}

/** Own scrolling, resizing and keyed DOM rows; callers supply row/spacer elements. */
function createVirtualRows({ viewport, container, key, createRow, updateRow,
  disposeRow = () => {}, createSpacer, createEmpty, headerHeight = () => 0, onRender = () => {} }) {
  const document = container.ownerDocument, view = document.defaultView;
  if (!viewport.hasAttribute("tabindex")) viewport.tabIndex = -1;
  const spacer = () => {
    const node = createSpacer ? createSpacer() : document.createElement("div");
    node.setAttribute("aria-hidden", "true"); container.append(node); return node;
  };
  const top = spacer(), bottom = spacer(), mounted = new Map();
  const empty = createEmpty?.();
  if (empty) container.insertBefore(empty, bottom);
  let rows = [], first = -1, last = -1, frame = null, destroyed = false, rendering = false, pending = false;
  const render = (force = false) => {
    if (destroyed) return;
    if (rendering) { pending = true; return; }
    const window = virtualRowWindow(rows.length, Math.max(0, viewport.scrollTop - headerHeight()), viewport.clientHeight);
    if (!force && first === window.first && last === window.last) return;
    const scrollTop = viewport.scrollTop;
    rendering = true;
    try {
      first = window.first; last = window.last;
      const visible = rows.slice(first, last), ids = new Set(visible.map(key));
      const remove = (node) => {
        const focused = container.getRootNode().activeElement || document.activeElement;
        if (node.contains(focused)) { focused.blur(); viewport.focus({ preventScroll: true }); }
        disposeRow(node); node.remove();
      };
      for (const [id, entry] of mounted) if (!ids.has(id)) { remove(entry.node); mounted.delete(id); }
      if (empty) empty.hidden = rows.length > 0;
      let next = empty || bottom;
      for (let index = visible.length - 1; index >= 0; index--) {
        const row = visible[index], id = key(row);
        let entry = mounted.get(id);
        if (!entry) { entry = { node: createRow(row) }; mounted.set(id, entry); }
        else updateRow?.(entry.node, row);
        entry.node.setAttribute("aria-rowindex", first + index + 2);
        if (entry.node.nextSibling !== next) container.insertBefore(entry.node, next);
        next = entry.node;
      }
      // A table spacer uses its cell for height; a list spacer uses the element.
      (top.firstElementChild || top).style.height = `${window.before}px`;
      (bottom.firstElementChild || bottom).style.height = `${window.after}px`;
      // Removing mounted rows can temporarily shrink the scroll area and clamp
      // its position, especially when jumping to the end of a native table.
      viewport.scrollTop = scrollTop;
      onRender([...mounted.values()].map((entry) => entry.node));
    } finally {
      rendering = false;
      if (pending && !destroyed) { pending = false; view.queueMicrotask(() => render(true)); }
    }
  };
  const schedule = () => { if (frame === null) frame = view.requestAnimationFrame(() => { frame = null; render(); }); };
  viewport.addEventListener("scroll", schedule);
  const observer = new view.ResizeObserver(schedule); observer.observe(viewport);
  return {
    setRows(next) { rows = next; render(true); },
    refresh() { render(true); },
    destroy() {
      destroyed = true; observer.disconnect(); viewport.removeEventListener("scroll", schedule);
      if (frame !== null) view.cancelAnimationFrame(frame);
      for (const { node } of mounted.values()) { disposeRow(node); node.remove(); }
      top.remove(); bottom.remove(); empty?.remove(); mounted.clear();
    },
  };
}

const DATA_EDITOR_CSS = `
${uiAppearanceStyles(":where(.cm-ui)")}
.cm-data-ui summary:hover { color:var(--cm-accent,var(--_cm-accent,#1677ff)); }
.cm-data-ui .cm-data-icon-button { display:inline-grid; place-items:center; width:var(--cm-icon-size,28px); min-width:var(--cm-icon-size,28px); min-height:var(--cm-icon-size,28px); padding:0; border:0; background:transparent; color:var(--cm-muted,var(--_cm-muted,#64748b)); line-height:1; }
.cm-data-ui .cm-data-icon-button:hover:not(:disabled) { background:var(--cm-hover,var(--_cm-hover,#e3ecf8)); color:var(--cm-accent,var(--_cm-accent,#1677ff)); }
.cm-data-icon-button svg { width:16px; height:16px; fill:none; stroke:currentColor; stroke-width:1.8; stroke-linecap:round; stroke-linejoin:round; }
.cm-data-toolbar { display:flex; flex-wrap:wrap; align-items:center; gap:8px; flex-shrink:0; }
.cm-data-error { color:var(--cm-danger,var(--_cm-danger,#b42318)); font-size:12px; margin:6px 0; }
.cm-data-error:empty { display:none; }
.cm-data-panel { display:flex; flex-direction:column; flex:1; min-width:0; min-height:0; gap:8px; }
.cm-data-tools { flex-shrink:0; }
.cm-data-actions { display:flex; align-items:center; gap:6px; flex-shrink:0; }
.cm-data-head { box-sizing:border-box; height:32px; border-bottom:1px solid var(--cm-border,var(--_cm-border,#dbe2eb)); background:var(--cm-subtle,var(--_cm-subtle,#f3f6fa)); color:var(--cm-muted,var(--_cm-muted,#556274)); font-size:11px; line-height:20px; font-weight:650; text-align:left; }
.cm-data-row { --cm-row-background:var(--cm-surface,var(--_cm-surface,#fff)); box-sizing:border-box; height:${DATA_ROW_HEIGHT}px; border-bottom:1px solid var(--cm-border,var(--_cm-border,#eee)); background:var(--cm-row-background); font-size:13px; line-height:20px; }
.cm-data-row > :is(th,td) { height:${DATA_ROW_HEIGHT}px; box-sizing:border-box; background:var(--cm-row-background); border-bottom:1px solid var(--cm-border,var(--_cm-border,#eee)); }
.cm-data-row:hover { --cm-row-background:var(--cm-hover,var(--_cm-hover,#e3ecf8)); }
.cm-data-row.is-selected { --cm-row-background:var(--cm-selection,var(--_cm-selection,#d6e6fc)); }
.cm-data-row:is(:hover,.is-selected) { --cm-edit-icon-width:20px; }
.cm-data-cell { box-sizing:border-box; padding:0 8px; text-align:left; vertical-align:middle; min-width:0; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; }
.cm-data-cell--select { width:32px; min-width:32px; padding:0; text-align:center; }
span.cm-data-cell--select { display:flex; align-items:center; justify-content:center; }
.cm-data-cell--select input { display:block; margin:auto; }
.cm-data-cell > input[type=color] { display:block; }
`;

function installDataEditorStyle(root) { installStyle(root, "cm-data-editor-style", (root.host ? controlStyles : lightControlStyles) + DATA_EDITOR_CSS); }

const [markup$a, css$a] = await Promise.all([
  "<details>\n  <summary part=\"trigger\"></summary>\n  <div part=\"panel\"><slot></slot></div>\n</details>\n",
  ":host { display: block; position: relative; color: var(--cm-foreground,var(--_cm-foreground, #273444)); font: inherit; }\n:host([hidden]) { display: none !important; }\nsummary {\n  display: block;\n  box-sizing: border-box;\n  min-height: var(--cm-control-height,var(--_cm-control-height,30px));\n  padding: var(--cm-control-padding,var(--_cm-control-padding,4px 8px));\n  border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de));\n  border-radius: var(--cm-control-radius,var(--_cm-control-radius,6px));\n  background: var(--cm-surface,var(--_cm-surface, #fff));\n  list-style: none;\n  white-space: nowrap;\n  cursor: pointer;\n}\nsummary::-webkit-details-marker { display: none; }\nsummary:hover, details[open] > summary { background: var(--cm-hover,var(--_cm-hover, #e3ecf8)); border-color: var(--cm-accent,var(--_cm-accent, #78a9ed)); }\nsummary:focus-visible { outline: 2px solid var(--cm-focus-ring,var(--_cm-focus-ring, #1677ff)); outline-offset: 2px; }\ndiv {\n  position: absolute;\n  right: 0;\n  top: calc(100% + 4px);\n  z-index: 10;\n  display: flex;\n  flex-direction: column;\n  gap: 8px;\n  width: max-content;\n  max-width: min(340px, 80vw);\n  max-height: 50vh;\n  overflow: auto;\n  padding: 10px;\n  border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de));\n  border-radius: 8px;\n  background: var(--cm-surface,var(--_cm-surface, #fff));\n  box-shadow: 0 4px 16px #0002;\n}\nslot { display: contents; }\n::slotted(button) { text-align: left; }\n::slotted(label) { display: flex; align-items: center; gap: 6px; }\n::slotted(input[type=number]) { width: 112px; }\n",
]);
const template$a = createTemplate(markup$a);
const activeMenus = new WeakMap();
const HTMLElementBase$a = globalThis.HTMLElement || class {};

class ActionMenu extends HTMLElementBase$a {
  static observedAttributes = ["label", "open"];
  #details;
  #trigger;
  #outside = (event) => { if (!event.composedPath().includes(this)) this.open = false; };

  constructor() {
    super();
    const root = mountComponent(this, template$a, css$a);
    this.#details = root.querySelector("details");
    this.#trigger = root.querySelector("summary");
    this.#details.ontoggle = () => { this.open = this.#details.open; };
    this.addEventListener("keydown", (event) => {
      if (event.key !== "Escape" || event.defaultPrevented || !this.open) return;
      this.open = false; this.#trigger.focus(); event.preventDefault(); event.stopPropagation();
    });
  }

  get open() { return this.#details.open; }
  set open(value) {
    this.#details.open = Boolean(value);
    this.toggleAttribute("open", this.#details.open);
  }
  connectedCallback() { this.#syncOpen(); }
  disconnectedCallback() { this.open = false; this.#release(); }

  attributeChangedCallback(name) {
    if (name === "label") this.#trigger.textContent = this.getAttribute("label") || "Actions";
    else { this.#details.open = this.hasAttribute("open"); this.#syncOpen(); }
  }

  #syncOpen() {
    if (!this.open || !this.isConnected) { this.#release(); return; }
    const previous = activeMenus.get(this.ownerDocument);
    if (previous && previous !== this) previous.open = false;
    activeMenus.set(this.ownerDocument, this);
    this.ownerDocument.addEventListener("pointerdown", this.#outside, true);
  }

  #release() {
    this.ownerDocument.removeEventListener("pointerdown", this.#outside, true);
    if (activeMenus.get(this.ownerDocument) === this) activeMenus.delete(this.ownerDocument);
  }
}

function defineActionMenu() {
  if (!customElements.get("clinker-action-menu")) customElements.define("clinker-action-menu", ActionMenu);
}

function selectionRange(orderedIds, anchor, id, range = false) {
  const start = orderedIds.indexOf(anchor), end = orderedIds.indexOf(id);
  const ranged = range && start >= 0 && end >= 0;
  return {
    ids: ranged ? orderedIds.slice(Math.min(start, end), Math.max(start, end) + 1) : [id],
    anchor: range && start >= 0 ? anchor : id,
    ranged,
  };
}

/** Selection shared by flat rows, tree leaves and annotation matrices. */
class RowSelection extends Set {
  #anchor = null;

  select(id, checked, { orderedIds = [], range = false, idsFor = (id) => [id] } = {}) {
    const { ids, anchor } = selectionRange(orderedIds, this.#anchor, id, range);
    this.setMany(ids.flatMap(idsFor), checked);
    this.#anchor = anchor;
  }

  setMany(ids, checked = true) {
    for (const id of ids) { if (checked) this.add(id); else this.delete(id); }
    this.#anchor = null;
  }

  prune(ids, anchorIds = ids) {
    const valid = ids instanceof Set ? ids : new Set(ids);
    for (const id of this) if (!valid.has(id)) this.delete(id);
    const anchors = anchorIds === ids ? valid : new Set(anchorIds);
    if (!anchors.has(this.#anchor)) this.#anchor = null;
  }

  clear() { super.clear(); this.#anchor = null; }
}

const [markup$9, css$9] = await Promise.all([
  "<span part=\"label\"></span>\n<button type=\"button\" part=\"edit-button\">\n  <svg aria-hidden=\"true\" viewBox=\"0 0 24 24\">\n    <path d=\"m4 20 4.1-1 10.2-10.2a2.1 2.1 0 0 0-3-3L5.1 16 4 20Z\"></path>\n    <path d=\"m13.8 7.2 3 3\"></path>\n  </svg>\n</button>\n<input type=\"text\" part=\"input\" hidden>\n<select part=\"input\" hidden></select>\n",
  ":host {\n  --_cm-control-height:28px;\n  display: flex;\n  align-items: center;\n  min-width: 0;\n  height: var(--cm-control-height,var(--_cm-control-height,28px));\n  box-sizing: border-box;\n  color: inherit;\n  font: inherit;\n}\n:host([hidden]), [hidden] { display: none !important; }\n:host([compact]) { flex: 0 0 auto; width: var(--cm-edit-width, max-content); max-width: 100%; }\n:host(:hover) { --cm-edit-icon-width:20px; }\nspan { flex: 0 1 auto; min-width: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; cursor: text; }\nspan[data-placeholder] { color: var(--cm-muted,var(--_cm-muted, #687586)); }\nbutton, input, select { box-sizing: border-box; font: inherit; color: inherit; }\nbutton {\n  display: inline-grid;\n  place-items: center;\n  flex: 0 0 auto;\n  width: var(--cm-edit-icon-width, 0px);\n  min-width: 0;\n  height: 20px;\n  min-height: 0;\n  overflow: hidden;\n  padding: 0;\n  border: 0;\n  border-radius: var(--cm-control-radius,var(--_cm-control-radius, 6px));\n  background: transparent;\n  color: var(--cm-muted,var(--_cm-muted, #64748b));\n  cursor: pointer;\n}\nbutton:hover { background: var(--cm-hover,var(--_cm-hover, #e3ecf8)); color: var(--cm-accent,var(--_cm-accent, #1677ff)); }\nbutton:focus-visible { width:20px; outline: 2px solid var(--cm-focus-ring,var(--_cm-focus-ring, #1677ff)); outline-offset: 2px; }\nbutton:disabled { cursor: default; opacity: .45; }\nsvg { width: 14px; height: 14px; fill: none; stroke: currentColor; stroke-width: 1.8; stroke-linecap: round; stroke-linejoin: round; }\ninput, select {\n  width: 100%;\n  min-width: 0;\n  height: var(--cm-control-height,var(--_cm-control-height,28px));\n  margin-inline: -1px;\n  padding: 2px 7px;\n  border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de));\n  border-radius: var(--cm-control-radius,var(--_cm-control-radius, 6px));\n  background: var(--cm-surface,var(--_cm-surface, #fff));\n  line-height: 20px;\n}\ninput::placeholder { color: var(--cm-muted,var(--_cm-muted, #687586)); }\ninput:focus, select:focus { outline: 2px solid var(--cm-focus-ring,var(--_cm-focus-ring, #1677ff)); outline-offset: -2px; }\n",
]);
const template$9 = createTemplate(markup$9);
const HTMLElementBase$9 = globalThis.HTMLElement || class {};

/** A field owns its draft; consumers accept or cancel value-commit events. */
class EditableField extends HTMLElementBase$9 {
  static observedAttributes = ["label", "placeholder", "type", "min", "max", "step", "required", "disabled"];
  #value = "";
  #revision = 0;
  #editing = false;
  #text;
  #button;
  #input;
  #select;
  #options = null;
  get #control() { return this.#options ? this.#select : this.#input; }

  constructor() {
    super();
    const root = mountComponent(this, template$9, css$9);
    this.#text = root.querySelector("span");
    this.#button = root.querySelector("button");
    this.#input = root.querySelector("input");
    this.#select = root.querySelector("select");
    this.#text.onclick = () => this.#beginEdit();
    this.#button.onclick = () => this.#beginEdit();
    for (const control of [this.#input, this.#select]) {
      control.onblur = () => { if (this.isConnected) this.#commit(); };
      control.oninput = () => control.setCustomValidity("");
      control.onkeydown = (event) => {
        if (event.isComposing) return;
        if (event.key === "Escape" || event.key === "Enter") {
          event.preventDefault(); event.stopPropagation();
          if (event.key === "Escape") this.#finish();
          else this.#commit();
          if (!this.#editing) this.#button.focus({ preventScroll: true });
        }
      };
    }
    this.#select.onchange = () => this.#commit();
    this.#render();
  }

  get value() { return this.#value; }
  set value(value) {
    const next = String(value ?? ""); this.#revision++;
    if (next === this.#value) return;
    this.#value = next;
    this.#render();
  }

  set options(entries) {
    this.#options = entries;
    if (entries) updateSelectOptions(this.#select, entries);
    this.#render();
  }

  attributeChangedCallback(name) {
    if (name === "disabled" && this.hasAttribute("disabled")) this.#finish();
    else this.#render();
  }
  disconnectedCallback() { this.#finish(); }

  cancelEdit() { this.#finish(); }
  setCustomValidity(message) { this.#control.setCustomValidity(message); }
  reportValidity() { return this.#control.reportValidity(); }

  #render() {
    const label = this.getAttribute("label") || "Value", placeholder = this.getAttribute("placeholder") || "";
    const display = this.#options?.find(([value]) => String(value) === this.#value)?.[1] ?? this.#value;
    this.#text.textContent = display || placeholder; this.#text.title = display || placeholder;
    this.#text.toggleAttribute("data-placeholder", !this.#value);
    this.#button.setAttribute("aria-label", `Edit ${label}`); this.#button.title = `Edit ${label}`;
    this.#select.setAttribute("aria-label", label);
    this.#input.setAttribute("aria-label", label); this.#input.placeholder = placeholder;
    const type = this.getAttribute("type") || "text";
    if (this.#input.type !== type) this.#input.type = type;
    for (const name of ["min", "max", "step"]) {
      if (this.hasAttribute(name)) this.#input.setAttribute(name, this.getAttribute(name));
      else this.#input.removeAttribute(name);
    }
    this.#select.required = this.hasAttribute("required");
    this.#input.required = this.hasAttribute("required");
    this.#select.disabled = this.#input.disabled = this.#button.disabled = this.hasAttribute("disabled");
    this.#text.hidden = this.#button.hidden = this.#editing;
    this.#input.hidden = !this.#editing || Boolean(this.#options);
    this.#select.hidden = !this.#editing || !this.#options;
    if (!this.#editing) this.#control.value = this.#value;
  }

  #beginEdit() {
    if (this.hasAttribute("disabled")) return;
    if (this.hasAttribute("compact")) this.style.setProperty("--cm-edit-width", `${this.getBoundingClientRect().width}px`);
    this.#editing = true;
    this.#control.value = this.#value; this.#control.setCustomValidity("");
    this.#render(); this.#control.focus({ preventScroll: true }); this.#control.select?.();
  }

  #finish() {
    this.#editing = false; this.#control.setCustomValidity("");
    this.style.removeProperty("--cm-edit-width"); this.#render();
  }

  #commit() {
    if (!this.#editing) return;
    if (!this.#control.reportValidity()) return;
    const value = this.#control.value, revision = this.#revision;
    if (value !== this.#value) {
      const event = new CustomEvent("value-commit", { detail: { value }, bubbles: true, composed: true, cancelable: true });
      if (!this.dispatchEvent(event)) { this.#control.reportValidity(); return; }
      if (revision === this.#revision) this.#value = value;
    }
    this.#finish();
  }
}

function defineEditableField() {
  if (!customElements.get("clinker-editable-field")) customElements.define("clinker-editable-field", EditableField);
}

const [markup$8, css$8] = await Promise.all([
  "<details>\n  <summary></summary>\n  <div data-items>\n    <template>\n      <div class=\"category\">\n        <input type=\"color\">\n        <clinker-editable-field required></clinker-editable-field>\n        <button type=\"button\" data-remove>×</button>\n      </div>\n    </template>\n    <form class=\"category\">\n      <input type=\"color\" aria-label=\"New category colour\">\n      <input type=\"text\" aria-label=\"New category name\" placeholder=\"New category\" required>\n      <button type=\"submit\" aria-label=\"Add category\">+</button>\n    </form>\n  </div>\n  <div data-removal hidden>\n    <span></span>\n    <select aria-label=\"Replacement category\"></select>\n    <button type=\"button\" data-confirm>Remove</button>\n    <button type=\"button\" data-cancel>Cancel</button>\n  </div>\n</details>\n",
  ":host { display: block; color: inherit; font: inherit; }\nsummary { cursor: pointer; color: var(--cm-muted,var(--_cm-muted, #687586)); }\nsummary:hover { color: var(--cm-accent,var(--_cm-accent, #1677ff)); }\n[data-items] { display: flex; flex-wrap: wrap; gap: 6px 16px; padding: 8px 0 2px; }\n.category, [data-removal] { display: flex; align-items: center; gap: 5px; }\n.category:hover { --cm-edit-icon-width:20px; }\nclinker-editable-field, input[type=text] { width: 140px; min-width: 0; }\n[data-removal] { flex-wrap: wrap; margin-top: 6px; padding: 6px; border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de)); border-radius: 4px; }\n",
]);
const template$8 = createTemplate(markup$8);
const HTMLElementBase$8 = globalThis.HTMLElement || class {};

class CategoryEditor extends HTMLElementBase$8 {
  #value = [];
  #palette = ["#888888"];
  #revision = 0;
  #rows = new Map();
  #details;
  #items;
  #template;
  #form;
  #name;
  #colour;
  #removal;
  #replacement;
  #removing = null;

  constructor() {
    super();
    const root = mountComponent(this, template$8, css$8);
    this.#details = root.querySelector("details"); this.#items = root.querySelector("[data-items]");
    this.#template = this.#items.querySelector("template"); this.#form = root.querySelector("form");
    this.#name = this.#form.querySelector('[type="text"]'); this.#colour = this.#form.querySelector('[type="color"]');
    this.#removal = root.querySelector("[data-removal]"); this.#replacement = this.#removal.querySelector("select");
    this.#name.oninput = () => this.#name.setCustomValidity("");
    this.#form.onsubmit = (event) => {
      event.preventDefault(); const name = this.#name.value.trim();
      if (this.#validateName(name, this.#name) && this.#propose([...this.#value, { name, colour: this.#colour.value }])) this.#name.value = "";
    };
    root.querySelector("[data-cancel]").onclick = () => this.#cancelRemoval();
    root.querySelector("[data-confirm]").onclick = () => {
      const replacement = this.#replacement.value === "missing" ? null : this.#value[Number(this.#replacement.value)].name;
      if (this.#propose(this.#value.filter(({ name }) => name !== this.#removing), [[this.#removing, replacement]])) {
        this.#removing = null; this.#renderRemoval(); this.#name.focus();
      }
    };
    this.addEventListener("keydown", (event) => {
      if (event.key === "Escape" && !event.defaultPrevented && this.#removing !== null) {
        event.preventDefault(); event.stopPropagation(); this.#cancelRemoval();
      }
    });
    this.#render();
  }

  get value() { return this.#value; }
  set value(items) {
    if (items.length === this.#value.length && items.every(({ name, colour }, index) =>
      name === this.#value[index].name && colour === this.#value[index].colour)) return;
    this.#value = items.map(({ name, colour }) => ({ name, colour }));
    this.#revision++; this.#render();
  }
  set palette(colours) { this.#palette = colours.length ? colours.map((colour) => colourHex(colour)) : ["#888888"]; this.#suggestColour(); }
  get open() { return this.#details.open; }
  set open(value) { this.#details.open = Boolean(value); }
  disconnectedCallback() { this.cancelEdit(); }
  cancelEdit() { for (const row of this.#rows.values()) row.querySelector("clinker-editable-field").cancelEdit(); }

  #validateName(name, control, previous) {
    control.setCustomValidity(!name ? "Category names cannot be empty." : name !== previous && this.#value.some((item) => item.name === name) ? "Category names must be unique." : "");
    return control.reportValidity();
  }

  #propose(value, replacements = []) {
    const revision = this.#revision;
    if (!this.dispatchEvent(new CustomEvent("value-commit", { detail: { value, replacements }, bubbles: true, composed: true, cancelable: true }))) return false;
    if (this.isConnected && revision === this.#revision) this.value = value;
    return true;
  }

  #createRow(name) {
    const row = this.ownerDocument.importNode(this.#template.content, true).firstElementChild;
    const colour = row.querySelector("input"), label = row.querySelector("clinker-editable-field"), remove = row.querySelector("button");
    colour.onchange = () => {
      if (!this.#propose(this.#value.map((item) => item.name === name ? { ...item, colour: colour.value } : item))) colour.value = colourHex(this.#value.find((item) => item.name === name).colour);
    };
    label.addEventListener("value-commit", (event) => {
      event.stopPropagation(); const next = event.detail.value.trim();
      if (!this.#validateName(next, label, name) || !this.#propose(this.#value.map((item) => item.name === name ? { ...item, name: next } : item), [[name, next]])) event.preventDefault();
      else label.value = next;
    });
    remove.onclick = () => { this.#removing = name; this.#renderRemoval(); this.#replacement.focus(); };
    colour.setAttribute("aria-label", `Colour for ${name}`); label.setAttribute("label", `Rename category ${name}`);
    remove.setAttribute("aria-label", `Remove category ${name}`);
    return row;
  }

  #render() {
    this.shadowRoot.querySelector("summary").textContent = `Categories (${this.#value.length})`;
    for (const [name, row] of this.#rows) if (!this.#value.some((item) => item.name === name)) {
      row.querySelector("clinker-editable-field").cancelEdit(); row.remove(); this.#rows.delete(name);
    }
    let next = this.#form;
    for (const item of [...this.#value].reverse()) {
      let row = this.#rows.get(item.name);
      if (!row) { row = this.#createRow(item.name); this.#rows.set(item.name, row); }
      row.querySelector("input").value = colourHex(item.colour); row.querySelector("clinker-editable-field").value = item.name;
      if (row.nextSibling !== next) this.#items.insertBefore(row, next);
      next = row;
    }
    if (!this.#value.some(({ name }) => name === this.#removing)) this.#removing = null;
    this.#renderRemoval(); this.#suggestColour();
  }

  #suggestColour() {
    const used = new Set(this.#value.map(({ colour }) => colourHex(colour)));
    this.#colour.value = this.#palette.find((colour) => !used.has(colour)) || this.#palette[this.#value.length % this.#palette.length] || "#888888";
  }

  #renderRemoval() {
    this.#removal.hidden = this.#removing === null;
    if (this.#removing === null) return;
    this.#removal.querySelector("span").textContent = `Replace “${this.#removing}” with`;
    updateSelectOptions(this.#replacement, [["missing", "Missing"], ...this.#value.flatMap(({ name }, index) => name === this.#removing ? [] : [[String(index), name]])]);
    this.#replacement.value = "missing";
  }

  #cancelRemoval() {
    const name = this.#removing; this.#removing = null; this.#renderRemoval();
    this.#rows.get(name)?.querySelector("button").focus();
  }
}

function defineCategoryEditor() {
  defineEditableField();
  if (!customElements.get("clinker-category-editor")) customElements.define("clinker-category-editor", CategoryEditor);
}

const [markup$7, css$7] = await Promise.all([
  "<select aria-label=\"Annotation panel\"></select>\n<button type=\"button\" data-action=\"add\">Add panel</button>\n<button type=\"button\" data-action=\"replace\">Replace data…</button>\n<button type=\"button\" data-action=\"show\">Hide</button>\n<button type=\"button\" data-action=\"left\" aria-label=\"Move left\" title=\"Move left\">←</button>\n<button type=\"button\" data-action=\"right\" aria-label=\"Move right\" title=\"Move right\">→</button>\n<button type=\"button\" data-action=\"remove\" aria-label=\"Remove panel\">Remove</button>\n<slot></slot>\n",
  ":host { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; flex-shrink: 0; color: inherit; font: inherit; }\n:host([embedded]), slot { display: contents; }\nselect { max-width: 260px; }\n:host([embedded]) select { max-width: 220px; }\n",
]);
const template$7 = createTemplate(markup$7);
const HTMLElementBase$7 = globalThis.HTMLElement || class {};

class AnnotationToolbar extends HTMLElementBase$7 {
  #select;
  #buttons;

  constructor() {
    super();
    const root = mountComponent(this, template$7, css$7);
    this.#select = root.querySelector("select");
    this.#buttons = Object.fromEntries([...root.querySelectorAll("button")].map((button) => [button.dataset.action, button]));
    this.#select.onchange = () => this.#emit("panel-change", { value: this.#select.value });
    root.addEventListener("click", (event) => {
      const button = event.target.closest("button[data-action]");
      if (button && !button.disabled) this.#emit("panel-action", { action: button.dataset.action });
    });
    this.state = {};
  }

  set state({ panels = [], value = null, creating = false }) {
    const panel = panels.find(({ uid }) => uid === value), index = panels.indexOf(panel);
    updateSelectOptions(this.#select, panels.map(({ uid, title }) => [uid, title || uid]));
    this.#select.value = value ?? ""; this.#select.hidden = !panels.length;
    for (const [action, button] of Object.entries(this.#buttons)) {
      button.hidden = action !== "add" && (!panel || creating || action === "replace" && panel.type === "groups");
    }
    this.#buttons.add.textContent = creating ? "Cancel" : "Add panel";
    this.#buttons.show.textContent = panel?.show === false ? "Show" : "Hide";
    this.#buttons.left.disabled = index <= 0;
    this.#buttons.right.disabled = index < 0 || index === panels.length - 1;
  }

  #emit(type, detail) { this.dispatchEvent(new CustomEvent(type, { detail, bubbles: true, composed: true })); }
}

function defineAnnotationToolbar() {
  if (!customElements.get("clinker-annotation-toolbar")) customElements.define("clinker-annotation-toolbar", AnnotationToolbar);
}

var EOL = {},
    EOF = {},
    QUOTE = 34,
    NEWLINE = 10,
    RETURN = 13;

function objectConverter(columns) {
  return new Function("d", "return {" + columns.map(function(name, i) {
    return JSON.stringify(name) + ": d[" + i + "] || \"\"";
  }).join(",") + "}");
}

function customConverter(columns, f) {
  var object = objectConverter(columns);
  return function(row, i) {
    return f(object(row), i, columns);
  };
}

// Compute unique columns in order of discovery.
function inferColumns(rows) {
  var columnSet = Object.create(null),
      columns = [];

  rows.forEach(function(row) {
    for (var column in row) {
      if (!(column in columnSet)) {
        columns.push(columnSet[column] = column);
      }
    }
  });

  return columns;
}

function pad(value, width) {
  var s = value + "", length = s.length;
  return length < width ? new Array(width - length + 1).join(0) + s : s;
}

function formatYear(year) {
  return year < 0 ? "-" + pad(-year, 6)
    : year > 9999 ? "+" + pad(year, 6)
    : pad(year, 4);
}

function formatDate(date) {
  var hours = date.getUTCHours(),
      minutes = date.getUTCMinutes(),
      seconds = date.getUTCSeconds(),
      milliseconds = date.getUTCMilliseconds();
  return isNaN(date) ? "Invalid Date"
      : formatYear(date.getUTCFullYear()) + "-" + pad(date.getUTCMonth() + 1, 2) + "-" + pad(date.getUTCDate(), 2)
      + (milliseconds ? "T" + pad(hours, 2) + ":" + pad(minutes, 2) + ":" + pad(seconds, 2) + "." + pad(milliseconds, 3) + "Z"
      : seconds ? "T" + pad(hours, 2) + ":" + pad(minutes, 2) + ":" + pad(seconds, 2) + "Z"
      : minutes || hours ? "T" + pad(hours, 2) + ":" + pad(minutes, 2) + "Z"
      : "");
}

function dsv(delimiter) {
  var reFormat = new RegExp("[\"" + delimiter + "\n\r]"),
      DELIMITER = delimiter.charCodeAt(0);

  function parse(text, f) {
    var convert, columns, rows = parseRows(text, function(row, i) {
      if (convert) return convert(row, i - 1);
      columns = row, convert = f ? customConverter(row, f) : objectConverter(row);
    });
    rows.columns = columns || [];
    return rows;
  }

  function parseRows(text, f) {
    var rows = [], // output rows
        N = text.length,
        I = 0, // current character index
        n = 0, // current line number
        t, // current token
        eof = N <= 0, // current token followed by EOF?
        eol = false; // current token followed by EOL?

    // Strip the trailing newline.
    if (text.charCodeAt(N - 1) === NEWLINE) --N;
    if (text.charCodeAt(N - 1) === RETURN) --N;

    function token() {
      if (eof) return EOF;
      if (eol) return eol = false, EOL;

      // Unescape quotes.
      var i, j = I, c;
      if (text.charCodeAt(j) === QUOTE) {
        while (I++ < N && text.charCodeAt(I) !== QUOTE || text.charCodeAt(++I) === QUOTE);
        if ((i = I) >= N) eof = true;
        else if ((c = text.charCodeAt(I++)) === NEWLINE) eol = true;
        else if (c === RETURN) { eol = true; if (text.charCodeAt(I) === NEWLINE) ++I; }
        return text.slice(j + 1, i - 1).replace(/""/g, "\"");
      }

      // Find next delimiter or newline.
      while (I < N) {
        if ((c = text.charCodeAt(i = I++)) === NEWLINE) eol = true;
        else if (c === RETURN) { eol = true; if (text.charCodeAt(I) === NEWLINE) ++I; }
        else if (c !== DELIMITER) continue;
        return text.slice(j, i);
      }

      // Return last token before EOF.
      return eof = true, text.slice(j, N);
    }

    while ((t = token()) !== EOF) {
      var row = [];
      while (t !== EOL && t !== EOF) row.push(t), t = token();
      if (f && (row = f(row, n++)) == null) continue;
      rows.push(row);
    }

    return rows;
  }

  function preformatBody(rows, columns) {
    return rows.map(function(row) {
      return columns.map(function(column) {
        return formatValue(row[column]);
      }).join(delimiter);
    });
  }

  function format(rows, columns) {
    if (columns == null) columns = inferColumns(rows);
    return [columns.map(formatValue).join(delimiter)].concat(preformatBody(rows, columns)).join("\n");
  }

  function formatBody(rows, columns) {
    if (columns == null) columns = inferColumns(rows);
    return preformatBody(rows, columns).join("\n");
  }

  function formatRows(rows) {
    return rows.map(formatRow).join("\n");
  }

  function formatRow(row) {
    return row.map(formatValue).join(delimiter);
  }

  function formatValue(value) {
    return value == null ? ""
        : value instanceof Date ? formatDate(value)
        : reFormat.test(value += "") ? "\"" + value.replace(/"/g, "\"\"") + "\""
        : value;
  }

  return {
    parse: parse,
    parseRows: parseRows,
    format: format,
    formatBody: formatBody,
    formatRows: formatRows,
    formatRow: formatRow,
    formatValue: formatValue
  };
}

var csv = dsv(",");
var csvParseRows = csv.parseRows;
var csvFormatRows = csv.formatRows;

var tsv = dsv("\t");
var tsvParseRows = tsv.parseRows;

const [markup$6, css$6] = await Promise.all([
  "<form>\n  <div data-create>\n    <label><span>Type</span><select name=\"type\" aria-label=\"Type\"><option value=\"heatmap\">Heatmap</option><option value=\"bar\">Bar chart</option><option value=\"groups\">Group markers</option></select></label>\n    <label><span>Title</span><input name=\"title\" type=\"text\" aria-label=\"Title\"></label>\n    <label data-scale><span>Values</span><select name=\"scale\" aria-label=\"Values\"><option value=\"continuous\">Numeric</option><option value=\"categorical\">Categories</option></select></label>\n    <p data-help></p>\n  </div>\n  <label class=\"data\"><span>Data</span><textarea name=\"data\" aria-label=\"Data\" rows=\"8\" spellcheck=\"false\"></textarea></label>\n  <p data-error role=\"alert\" hidden></p>\n  <div class=\"actions\"><button type=\"submit\">Add panel</button><button type=\"button\" data-cancel hidden>Cancel</button></div>\n</form>\n",
  ":host { display: block; flex-shrink: 0; min-width: 0; color: inherit; font: inherit; }\nlabel { display: flex; align-items: center; gap: 8px; margin-bottom: 8px; }\nlabel > span { flex: 1; min-width: 0; }\ninput, select { width: 140px; max-width: 65%; }\n.data { display: block; }\ntextarea { display: block; width: 100%; margin-top: 4px; resize: vertical; }\np { margin: 8px 0; color: var(--cm-muted,var(--_cm-muted, #687586)); }\n[data-error] { color: var(--cm-danger,var(--_cm-danger, #b42318)); }\n.actions { display: flex; gap: 8px; }\n",
]);
const template$6 = createTemplate(markup$6);
const HTMLElementBase$6 = globalThis.HTMLElement || class {};

class AnnotationDataForm extends HTMLElementBase$6 {
  static observedAttributes = ["mode"];
  #form;
  #clusters = [];

  constructor() {
    super();
    const root = mountComponent(this, template$6, css$6);
    this.#form = root.querySelector("form");
    this.#field("type").onchange = () => { this.#sync(); this.#example(); };
    this.#form.onsubmit = (event) => { event.preventDefault(); this.error = ""; this.#emit("data-submit", this.value); };
    root.querySelector("[data-cancel]").onclick = () => this.#emit("data-cancel");
    this.#sync();
  }

  set clusters(value) { this.#clusters = value; this.#example(); }
  get value() { return Object.fromEntries(["type", "title", "scale", "data"].map((name) => [name, this.#field(name).value])); }
  set value(value) { for (const [name, text] of Object.entries(value)) this.#field(name).value = text; this.#sync(); }
  set error(message) {
    const error = this.shadowRoot.querySelector("[data-error]"); error.textContent = message; error.hidden = !message;
  }
  focus(options) { this.#field(this.getAttribute("mode") === "replace" ? "data" : "type").focus(options); }
  attributeChangedCallback() { this.#sync(); }

  #field(name) { return this.#form.elements.namedItem(name); }
  #emit(type, detail) { this.dispatchEvent(new CustomEvent(type, { detail, bubbles: true, composed: true })); }
  #sync() {
    const root = this.shadowRoot, replace = this.getAttribute("mode") === "replace", type = this.#field("type").value;
    root.querySelector("[data-create]").hidden = replace;
    root.querySelector("[data-scale]").hidden = type === "groups";
    this.#field("scale").disabled = type === "bar";
    if (type === "bar") this.#field("scale").value = "continuous";
    root.querySelector("[data-help]").textContent = type === "groups"
      ? "Define groups as JSON. Each needs a uid, a label, and either clusterIds or a tree nodeUid. Colours are optional."
      : "Paste CSV or a tab-separated table. First column: cluster ID; remaining headings: data columns. Blank cells are missing. Bars use one data column. JSON {columns, values} is also accepted.";
    this.#field("data").rows = replace ? 5 : 8;
    this.#field("data").setAttribute("aria-label", replace ? "Replacement annotation data (CSV, TSV or JSON)" : "Data");
    root.querySelector('[type="submit"]').textContent = replace ? "Replace data" : "Add panel";
    root.querySelector("[data-cancel]").hidden = !replace;
  }
  #example() {
    if (this.getAttribute("mode") === "replace") return;
    this.#field("data").value = this.#field("type").value === "groups"
      ? JSON.stringify({ groups: [{ uid: "group-1", label: "Example group", clusterIds: this.#clusters.slice(0, 2).map(({ uid }) => uid) }] }, null, 2)
      : csvFormatRows([["cluster", "Value"], ...this.#clusters.map(({ uid }) => [uid, ""])]);
  }
}

function defineAnnotationDataForm() {
  if (!customElements.get("clinker-annotation-data-form")) customElements.define("clinker-annotation-data-form", AnnotationDataForm);
}

const [markup$5, css$5] = await Promise.all([
  "<div data-groups>\n  <template>\n    <fieldset>\n      <legend></legend>\n      <div class=\"controls\">\n        <clinker-editable-field label=\"Group label\"></clinker-editable-field>\n        <input type=\"color\" aria-label=\"Group colour\">\n        <select data-mode aria-label=\"Membership\"></select>\n        <button type=\"button\">Remove group</button>\n      </div>\n      <select data-members></select>\n    </fieldset>\n  </template>\n</div>\n<button type=\"button\" data-add>Add group</button>\n",
  ":host { display: block; color: inherit; font: inherit; }\nfieldset { min-width: 0; margin: 0 0 8px; padding: 8px; border: 1px solid var(--cm-border,var(--_cm-border, #cbd3de)); border-radius: 6px; }\nlegend { padding: 0 4px; font-weight: 600; }\n.controls { display: flex; flex-wrap: wrap; align-items: center; gap: 8px; }\nclinker-editable-field { width: 140px; }\nclinker-editable-field:hover { --cm-edit-icon-width: 20px; }\nselect { max-width: 100%; }\n[data-members] { display: block; margin-top: 8px; }\n",
]);
const template$5 = createTemplate(markup$5);
const HTMLElementBase$5 = globalThis.HTMLElement || class {};

class AnnotationGroupEditor extends HTMLElementBase$5 {
  #state = { groups: [], clusters: [], tree: null };
  #revision = 0;
  #rows = new Map();
  #items;
  #template;
  #add;

  constructor() {
    super();
    const root = mountComponent(this, template$5, css$5);
    this.#items = root.querySelector("[data-groups]"); this.#template = this.#items.querySelector("template");
    this.#add = root.querySelector("[data-add]");
    this.#add.onclick = () => {
      const { groups, clusters } = this.#state;
      if (!clusters.length) return;
      let serial = 1; while (groups.some(({ uid }) => uid === `group-${serial}`)) serial++;
      this.#propose([...groups, { uid: `group-${serial}`, label: `Group ${serial}`, clusterIds: [clusters[0].uid] }]);
    };
    this.#render();
  }

  set state(value) { this.#state = value; this.#revision++; this.#render(); }
  cancelEdit() { for (const { label } of this.#rows.values()) label.cancelEdit(); }
  disconnectedCallback() { this.cancelEdit(); }

  #propose(groups) {
    const revision = this.#revision;
    const accepted = this.dispatchEvent(new CustomEvent("value-commit", { detail: { value: groups }, bubbles: true, composed: true, cancelable: true }));
    if (accepted && revision === this.#revision) this.state = { ...this.#state, groups };
    if (!accepted) this.#render();
    return accepted;
  }

  #edit(uid, changes) { return this.#propose(this.#state.groups.map((group) => group.uid === uid ? { ...group, ...changes } : group)); }

  #createRow(uid) {
    const node = this.ownerDocument.importNode(this.#template.content, true).firstElementChild;
    const row = { node, label: node.querySelector("clinker-editable-field"), colour: node.querySelector("input"),
      mode: node.querySelector("[data-mode]"), members: node.querySelector("[data-members]"), remove: node.querySelector("button") };
    row.label.addEventListener("value-commit", (event) => {
      event.stopPropagation(); if (!this.#edit(uid, { label: event.detail.value })) event.preventDefault();
    });
    row.colour.onchange = () => this.#edit(uid, { colour: row.colour.value });
    row.remove.onclick = () => this.#propose(this.#state.groups.filter((group) => group.uid !== uid));
    row.mode.onchange = () => {
      const { tree, clusters } = this.#state;
      if (row.mode.value === "tree" && tree) this.#edit(uid, { clusterIds: undefined, nodeUid: tree.root });
      else if (clusters.length) this.#edit(uid, { nodeUid: undefined, clusterIds: [clusters[0].uid] });
      else this.#render();
    };
    row.members.onchange = () => {
      const { clusters, tree } = this.#state;
      this.#edit(uid, row.members.multiple
        ? { clusterIds: [...row.members.selectedOptions].map((option) => clusters.find(({ uid }) => String(uid) === option.value).uid) }
        : { nodeUid: tree.nodes.find(({ uid }) => String(uid) === row.members.value).uid });
    };
    return row;
  }

  #render() {
    const { groups, clusters, tree, colour } = this.#state;
    this.#add.disabled = !clusters.length;
    for (const [uid, row] of this.#rows) if (!groups.some((group) => group.uid === uid)) {
      row.label.cancelEdit(); row.node.remove(); this.#rows.delete(uid);
    }
    let next = null;
    for (const group of [...groups].reverse()) {
      let row = this.#rows.get(group.uid);
      if (!row) { row = this.#createRow(group.uid); this.#rows.set(group.uid, row); }
      if (row.node.nextSibling !== next || row.node.parentNode !== this.#items) this.#items.insertBefore(row.node, next);
      next = row.node;
      row.node.querySelector("legend").textContent = group.label || group.uid;
      row.label.value = group.label ?? group.uid; row.colour.value = colourHex(group.colour ?? colour ?? "#267c9e");
      row.remove.disabled = groups.length === 1;
      updateSelectOptions(row.mode, [["clusters", "Selected clusters"], ...(tree ? [["tree", "Tree clade"]] : [])]);
      row.mode.value = group.nodeUid === undefined ? "clusters" : "tree";
      row.members.multiple = group.nodeUid === undefined; row.members.size = row.members.multiple ? Math.max(1, Math.min(5, clusters.length)) : 0;
      row.members.setAttribute("aria-label", row.members.multiple ? `Clusters in ${group.label || group.uid}` : "Tree node");
      updateSelectOptions(row.members, row.members.multiple ? clusters.map(({ uid, name }) => [uid, `${name || uid} (${uid})`]) : (tree?.nodes || []).map(({ uid, label }) => [uid, label || uid]));
      for (const option of row.members.options) option.selected = row.members.multiple
        ? (group.clusterIds || []).some((uid) => String(uid) === option.value) : String(group.nodeUid) === option.value;
    }
  }
}

function defineAnnotationGroupEditor() {
  defineEditableField();
  if (!customElements.get("clinker-annotation-group-editor")) customElements.define("clinker-annotation-group-editor", AnnotationGroupEditor);
}

const [markup$4, css$4] = await Promise.all([
  "<span part=\"summary\"></span>\n<div part=\"actions\" hidden>\n  <slot name=\"actions\"></slot>\n  <button type=\"button\" data-clear aria-label=\"Clear selection\">Clear</button>\n</div>\n<button type=\"button\" data-all hidden>Select all visible</button>\n",
  ":host {\n  display: flex;\n  align-items: center;\n  gap: 8px;\n  height: 32px;\n  min-width: 0;\n  flex-shrink: 0;\n  color: var(--cm-muted,var(--_cm-muted, #687586));\n  font: 12px/1.4 system-ui, sans-serif;\n}\nspan { flex: 1; min-width: 0; overflow: hidden; text-overflow: ellipsis; white-space: nowrap; }\nspan.is-selected { color: var(--cm-foreground,var(--_cm-foreground, #273444)); font-weight: 600; }\ndiv { display: flex; align-items: center; gap: 6px; flex-shrink: 0; }\nslot { display: contents; }\nbutton { flex-shrink: 0; }\n",
]);
const template$4 = createTemplate(markup$4);
const HTMLElementBase$4 = globalThis.HTMLElement || class {};

class SelectionBar extends HTMLElementBase$4 {
  #summary;
  #actions;
  #all;
  #ids = [];

  constructor() {
    super();
    const root = mountComponent(this, template$4, css$4);
    this.#summary = root.querySelector("span");
    this.#actions = root.querySelector("div");
    this.#all = root.querySelector("[data-all]");
    for (const [selector, name] of [["[data-clear]", "clear-selection"], ["[data-all]", "select-all"]]) {
      root.querySelector(selector).onclick = () => this.dispatchEvent(new CustomEvent(name, { detail: { ids: this.#ids }, bubbles: true, composed: true }));
    }
  }

  set state({ selection, ids, summary, noun, selectSuffix = "" }) {
    this.#ids = ids;
    const count = selection.size, total = ids.length;
    const all = total > 0 && ids.every((id) => selection.has(id));
    this.#summary.textContent = count ? `${count.toLocaleString()} ${noun}${count === 1 ? "" : "s"} selected` : summary;
    this.#summary.classList.toggle("is-selected", count > 0);
    this.#summary.title = summary;
    this.#actions.hidden = !count;
    this.#all.hidden = !total;
    this.#all.disabled = all;
    this.#all.textContent = `Select all ${total.toLocaleString()} visible${selectSuffix}`;
  }
}

function defineSelectionBar() {
  if (!customElements.get("clinker-selection-bar")) customElements.define("clinker-selection-bar", SelectionBar);
}

const [markup$3, css$3] = await Promise.all([
  "<div part=\"viewport\" tabindex=\"-1\"><slot name=\"header\"></slot><slot></slot><p part=\"empty\" hidden>No records match the filter.</p></div>\n<template id=\"heading\"><span class=\"cm-data-cell\"></span></template>\n<template id=\"sort\"><button type=\"button\" class=\"cm-data-cell cm-editor__sort\"><span data-label></span><span aria-hidden=\"true\"></span></button></template>\n",
  ":host { display: flex; flex: 1; min-width: 0; min-height: 120px; color: inherit; font: inherit; }\n:host([hidden]) { display: none !important; }\n[hidden] { display: none !important; }\n[part=viewport] { flex: 1; min-width: 0; min-height: 0; overflow: auto; overflow-anchor: none; outline: none; }\nslot { display: contents; }\n[part=empty] { padding: 16px; color: var(--cm-muted,var(--_cm-muted, #667386)); text-align: center; }\n",
]);
const template$3 = createTemplate(markup$3);
const HTMLElementBase$3 = globalThis.HTMLElement || class {};

class DataTable extends HTMLElementBase$3 {
  #viewport;
  #renderer = null;
  #rows = [];
  #view = null;
  #header = null;
  #columnsKey;

  constructor() {
    super();
    const root = mountComponent(this, template$3, css$3);
    this.#viewport = root.querySelector("div");
    this.addEventListener("mousedown", (event) => {
      if (event.button === 0 && event.shiftKey && event.target.closest(".cm-data-row") &&
          !event.target.closest('input:not([type="checkbox"]), textarea, select, [contenteditable], clinker-editable-field')) event.preventDefault();
    });
  }

  connectedCallback() { this.#mount(); }
  disconnectedCallback() { this.#dispose(); }

  // Row adapters provide their container, key and create/update/dispose callbacks.
  set renderer(value) {
    if (value === this.#renderer) return;
    this.#dispose(); this.#renderer = value; this.#rows = []; this.#syncEmpty(); this.#mount();
  }
  set rows(value) { this.#rows = value; this.#syncEmpty(); this.#view?.setRows(value); }
  set header({ columns, sort }) {
    const node = this.querySelector('[slot="header"]');
    if (!node) return;
    const key = JSON.stringify(columns);
    if (this.#header !== node || this.#columnsKey !== key) {
      this.#header = node; this.#columnsKey = key;
      node.replaceChildren(...columns.map(({ label, sort }, index) => {
        const template = this.shadowRoot.getElementById(sort ? "sort" : "heading");
        const cell = this.ownerDocument.importNode(template.content, true).firstElementChild;
        if (sort) { cell.dataset.sort = sort; cell.querySelector("[data-label]").textContent = label; cell.title = `Sort by ${label}`; }
        else cell.textContent = label;
        cell.classList.toggle("cm-data-cell--select", index === 0); return cell;
      }));
    }
    for (const button of node.querySelectorAll("[data-sort]")) {
      const label = button.querySelector("[data-label]").textContent, active = sort?.key === button.dataset.sort;
      button.setAttribute("aria-label", `Sort by ${label}, ${active && sort.direction < 0 ? "descending" : "ascending"}`);
      button.lastElementChild.textContent = active ? sort.direction < 0 ? "↓" : "↑" : "↕";
    }
  }
  get scrollTop() { return this.#viewport.scrollTop; }
  set scrollTop(value) { this.#viewport.scrollTop = value; }
  get scrollHeight() { return this.#viewport.scrollHeight; }
  scrollTo(...args) { this.#viewport.scrollTo(...args); }

  #mount() {
    if (!this.isConnected || !this.#renderer || this.#view) return;
    this.#view = createVirtualRows({ ...this.#renderer, viewport: this.#viewport });
    this.#view.setRows(this.#rows);
  }
  #dispose() { this.#view?.destroy(); this.#view = null; }
  #syncEmpty() { this.shadowRoot.querySelector("[part=empty]").hidden = !this.#renderer || this.#rows.length > 0; }
}

function defineDataTable() {
  if (!customElements.get("clinker-data-table")) customElements.define("clinker-data-table", DataTable);
}

const [markup$2, css$2] = await Promise.all([
  "<select data-columns aria-label=\"Columns to edit\"></select>\n<input type=\"number\" step=\"any\" placeholder=\"Value\" aria-label=\"Batch value\">\n<select data-categories aria-label=\"Batch value\" hidden></select>\n<button type=\"button\" data-set>Set values</button>\n<button type=\"button\" data-clear>Clear values</button>\n",
  ":host { display: flex; flex-direction: column; gap: 8px; color: inherit; font: inherit; }\nbutton { text-align: left; }\ninput { width: 112px; }\n",
]);
const template$2 = createTemplate(markup$2);
const HTMLElementBase$2 = globalThis.HTMLElement || class {};

class AnnotationValueActions extends HTMLElementBase$2 {
  #uid;
  #columns;
  #numeric;
  #categories;
  #domain = null;

  constructor() {
    super();
    const root = mountComponent(this, template$2, css$2);
    this.#columns = root.querySelector("[data-columns]"); this.#numeric = root.querySelector("input"); this.#categories = root.querySelector("[data-categories]");
    root.querySelector("[data-set]").onclick = () => this.#submit(false);
    root.querySelector("[data-clear]").onclick = () => this.#submit(true);
  }

  set state({ uid, columns, domain = null }) {
    const changed = this.#uid !== uid, value = changed ? "" : this.#input.value;
    this.#uid = uid; this.#domain = domain;
    this.#numeric.hidden = domain !== null; this.#categories.hidden = domain === null;
    updateSelectOptions(this.#columns, [["all", "All columns"], ...columns.map((name, index) => [String(index), name])]);
    if (domain !== null) updateSelectOptions(this.#categories, [["", "Missing"], ...domain.map((name, index) => [String(index), name])]);
    if (changed || !this.#columns.value) this.#columns.value = "0";
    this.#input.value = value;
    if (domain !== null && this.#categories.selectedIndex < 0) this.#categories.value = "";
  }

  get #input() { return this.#domain === null ? this.#numeric : this.#categories; }
  #submit(clear) {
    const input = this.#input;
    if (!clear && !input.reportValidity()) return;
    const value = clear || input.value === "" ? null : this.#domain === null ? Number(input.value) : this.#domain[Number(input.value)];
    const column = this.#columns.value === "all" ? "all" : Number(this.#columns.value);
    this.dispatchEvent(new CustomEvent("values-commit", { detail: { column, value }, bubbles: true, composed: true, cancelable: true }));
  }
}

function defineAnnotationValueActions() {
  if (!customElements.get("clinker-annotation-value-actions")) customElements.define("clinker-annotation-value-actions", AnnotationValueActions);
}

const [markup$1, css$1] = await Promise.all([
  "<div class=\"cm-ui cm-data-panel cm-data-ui\">\n  <slot name=\"tools\"></slot>\n  <clinker-selection-bar class=\"cm-data-meta\">\n    <div slot=\"actions\" class=\"cm-data-actions\">\n      <clinker-action-menu label=\"Edit values…\"><clinker-annotation-value-actions></clinker-annotation-value-actions></clinker-action-menu>\n    </div>\n  </clinker-selection-bar>\n  <clinker-data-table class=\"cm-data-table cm-annotation-matrix\" role=\"region\" aria-label=\"Editable annotation data\"></clinker-data-table>\n</div>\n<template id=\"table\">\n<table><thead><tr>\n  <th class=\"cm-data-cell cm-data-cell--select cm-data-head\"><input type=\"checkbox\" aria-label=\"Select all filtered clusters\"></th>\n  <th class=\"cm-data-cell cm-data-head\">Cluster</th>\n</tr></thead><tbody></tbody></table>\n</template>\n\n<template id=\"column\">\n<th class=\"cm-data-cell cm-data-head\"><div class=\"cm-annotation-column\">\n  <clinker-editable-field required></clinker-editable-field><button type=\"button\" class=\"cm-data-icon-button\">×</button>\n</div></th>\n</template>\n\n<template id=\"add-column\">\n<th class=\"cm-data-cell cm-data-head cm-annotation-add-column\"><form>\n  <button type=\"button\" aria-label=\"Add column\">+</button>\n  <input type=\"text\" aria-label=\"New column name\" placeholder=\"Column name\" hidden><button type=\"submit\" hidden>Add</button>\n</form></th>\n</template>\n\n<template id=\"row\">\n<tr class=\"cm-data-row\">\n  <td class=\"cm-data-cell cm-data-cell--select\"><input type=\"checkbox\"></td><th class=\"cm-data-cell\" scope=\"row\"></th>\n</tr>\n</template>\n\n<template id=\"empty-cell\">\n<td class=\"cm-data-cell\"></td>\n</template>\n\n<template id=\"spacer\">\n<tr class=\"cm-annotation-spacer\" aria-hidden=\"true\"><td></td></tr>\n</template>\n\n<template id=\"cell\">\n  <td class=\"cm-data-cell\"><clinker-editable-field placeholder=\"—\" step=\"any\"></clinker-editable-field></td>\n</template>\n",
  ":host { display:flex; flex:1; min-width:0; min-height:0; color:inherit; font:inherit; }\nslot { display:contents; }\n.cm-annotation-matrix { position:relative; }\n.cm-annotation-matrix table { border-collapse:separate; border-spacing:0; font:inherit; color:inherit; width:max-content; min-width:100%; }\n\n.cm-annotation-matrix thead th { position:sticky; top:0; z-index:2; background:var(--cm-subtle,var(--_cm-subtle,#f3f6fa)); }\n.cm-annotation-matrix tr > :first-child { position:sticky; left:0; z-index:1; width:32px; min-width:32px; background:var(--cm-row-background,var(--cm-subtle,#f3f6fa)); }\n.cm-annotation-matrix tr > :nth-child(2) { position:sticky; left:32px; z-index:1; }\n.cm-annotation-matrix tbody th { background:var(--cm-row-background); }\n.cm-annotation-matrix thead th:nth-child(-n+2) { z-index:3; background:var(--cm-subtle,var(--_cm-subtle,#f3f6fa)); }\n.cm-annotation-matrix tbody th { cursor:pointer; min-width:200px; max-width:340px; overflow:hidden; text-overflow:ellipsis; white-space:nowrap; font-weight:400; }\n.cm-annotation-matrix :is(input:not([type=checkbox]),clinker-editable-field) { width:140px; }\n.cm-annotation-column { display:flex; align-items:center; }\n.cm-annotation-column button { --cm-icon-size:24px; visibility:hidden; }\n.cm-annotation-column:is(:hover,:focus-within) button { visibility:visible; }\n.cm-annotation-add-column { width:36px; }\n.cm-annotation-add-column form { display:flex; gap:4px; }\n.cm-annotation-spacer > td { height:0; border:0; padding:0; }\n",
]);
const template$1 = createTemplate(markup$1);
const HTMLElementBase$1 = globalThis.HTMLElement || class {};

class AnnotationTable extends HTMLElementBase$1 {
  #viewport;
  #panel;
  #domain = null;
  #filtered = [];
  #filter = "";
  #schema;
  #table;
  #selection = new RowSelection();
  #rowIds = new WeakMap();
  #selectionBar;
  #batch;

  constructor() {
    super();
    const root = mountComponent(this, template$1, DATA_EDITOR_CSS + css$1);
    this.#viewport = root.querySelector("clinker-data-table");
    this.#selectionBar = root.querySelector("clinker-selection-bar");
    this.#selectionBar.addEventListener("clear-selection", () => this.clearSelection());
    this.#selectionBar.addEventListener("select-all", (event) => { this.#selection.setMany(event.detail.ids); this.#syncSelection(); });
    this.#batch = root.querySelector("clinker-annotation-value-actions");
    this.#batch.addEventListener("values-commit", (event) => {
      event.stopPropagation();
      const { column, value } = event.detail;
      if (!this.#setValues(this.#selection, column, value)) event.preventDefault();
    });
    root.addEventListener("keydown", (event) => {
      if (event.key === "Escape" && !event.defaultPrevented && this.#selection.size) {
        event.preventDefault(); event.stopPropagation(); this.clearSelection();
      }
    });
  }

  set state({ panel, clusters, domain = null, filter = "" }) {
    const changed = this.#panel?.uid !== panel.uid, filtered = this.#filter !== filter;
    if (changed) this.#selection.clear();
    this.#panel = panel; this.#domain = domain; this.#filter = filter;
    this.#selection.prune(new Set(clusters.map(({ uid }) => uid)));
    const scroll = changed || filtered ? 0 : this.#viewport.scrollTop;
    const schema = JSON.stringify([panel.uid, panel.type, panel.columns, domain]);
    if (this.#schema !== schema) { this.#schema = schema; this.#build(); }
    const needle = filter.trim().toLowerCase();
    this.#filtered = clusters.filter(({ uid, name }) => `${name} ${uid}`.toLowerCase().includes(needle));
    this.#viewport.rows = this.#filtered;
    this.#viewport.scrollTop = scroll;
  }

  clearSelection() { this.#selection.clear(); this.#syncSelection(); }
  cancelEdit() { for (const field of this.shadowRoot.querySelectorAll("clinker-editable-field")) field.cancelEdit(); }
  disconnectedCallback() { this.cancelEdit(); }

  #element(name) { return this.ownerDocument.importNode(this.shadowRoot.getElementById(name).content, true).firstElementChild; }
  #emit(type, detail) {
    const accepted = this.dispatchEvent(new CustomEvent(type, { detail, bubbles: true, composed: true, cancelable: true }));
    if (accepted) for (const menu of this.shadowRoot.querySelectorAll("clinker-action-menu")) menu.open = false;
    return accepted;
  }
  #columns(columns, values) { return this.#emit("columns-change", { columns, ...(values ? { values } : {}) }); }
  #setValues(ids, column, value) {
    const values = {};
    for (const uid of ids) {
      const row = this.#panel.values[uid] || this.#panel.columns.map(() => null);
      values[uid] = row.map((current, index) => column === "all" || column === index ? value : current);
    }
    return this.#emit("rows-change", { values });
  }

  #syncSelection() {
    this.#selectionBar.state = { selection: this.#selection, ids: this.#filtered.map(({ uid }) => uid), summary: `${this.#filtered.length.toLocaleString()} visible clusters`, noun: "cluster" };
    if (!this.#table) return;
    const all = this.#table.querySelector("thead input"), selected = this.#filtered.filter(({ uid }) => this.#selection.has(uid)).length;
    all.checked = this.#filtered.length > 0 && selected === this.#filtered.length; all.indeterminate = selected > 0 && !all.checked;
    for (const row of this.#table.querySelectorAll("tbody .cm-data-row")) {
      const checked = this.#selection.has(this.#rowIds.get(row));
      row.classList.toggle("is-selected", checked); row.setAttribute("aria-selected", String(checked)); row.querySelector("input").checked = checked;
    }
  }

  #build() {
    const viewport = this.#viewport, panel = this.#panel, extra = panel.type === "heatmap";
    this.cancelEdit(); viewport.renderer = null;
    this.#table = this.#element("table"); viewport.replaceChildren(this.#table);
    const head = this.#table.querySelector("tr");
    head.querySelector("input").onchange = (event) => { this.#selection.setMany(this.#filtered.map(({ uid }) => uid), event.target.checked); this.#syncSelection(); };
    this.#batch.state = { uid: panel.uid, columns: panel.columns, domain: this.#domain };
    panel.columns.forEach((column, index) => {
      const heading = this.#element("column"), field = heading.querySelector("clinker-editable-field"), remove = heading.querySelector("button");
      field.value = column; field.setAttribute("label", `Column ${index + 1} name`);
      field.addEventListener("value-commit", (event) => {
        event.stopPropagation(); const columns = [...this.#panel.columns]; columns[index] = event.detail.value.trim();
        if (!this.#columns(columns)) event.preventDefault();
      });
      if (extra) {
        remove.setAttribute("aria-label", `Remove column ${column}`); remove.disabled = panel.columns.length === 1;
        remove.onclick = () => this.#columns(this.#panel.columns.filter((_, i) => i !== index),
          Object.fromEntries(Object.entries(this.#panel.values).map(([uid, values]) => [uid, values.filter((_, i) => i !== index)])));
      } else remove.remove();
      head.append(heading);
    });
    if (extra) {
      const heading = this.#element("add-column"), form = heading.querySelector("form"), [add, save] = form.querySelectorAll("button"), name = form.querySelector("input");
      add.onclick = () => { add.hidden = true; name.hidden = save.hidden = false; name.focus(); };
      form.onsubmit = (event) => {
        event.preventDefault(); if (!name.value.trim()) return;
        this.#columns([...this.#panel.columns, name.value.trim()], Object.fromEntries(Object.entries(this.#panel.values).map(([uid, values]) => [uid, [...values, null]])));
      };
      name.onkeydown = (event) => { if (event.key === "Escape") { event.preventDefault(); name.hidden = save.hidden = true; add.hidden = false; add.focus(); } };
      head.append(heading);
    }
    viewport.renderer = { container: this.#table.querySelector("tbody"), key: ({ uid }) => uid,
      createRow: (cluster) => this.#createRow(cluster, extra), updateRow: (row, cluster) => this.#updateRow(row, cluster),
      disposeRow: (row) => { for (const field of row.querySelectorAll("clinker-editable-field")) field.cancelEdit(); },
      createSpacer: () => { const row = this.#element("spacer"); row.firstElementChild.colSpan = panel.columns.length + 2 + Number(extra); return row; },
      headerHeight: () => head.offsetHeight, onRender: () => { this.#table.setAttribute("aria-rowcount", this.#filtered.length + 1); this.#syncSelection(); } };
  }

  #createRow(cluster, extra) {
    const row = this.#element("row"), check = row.querySelector("input"), name = row.querySelector("th"); this.#rowIds.set(row, cluster.uid);
    const select = (checked, range) => { this.#selection.select(cluster.uid, checked, { range, orderedIds: this.#filtered.map(({ uid }) => uid) }); this.#syncSelection(); };
    check.onclick = (event) => select(check.checked, event.shiftKey);
    name.onclick = (event) => select(!this.#selection.has(cluster.uid), event.shiftKey);
    this.#panel.columns.forEach((_, index) => {
      const cell = this.#element("cell"), field = cell.firstElementChild;
      if (this.#domain === null) field.setAttribute("type", "number");
      else field.options = [["", "Missing"], ...this.#domain.map((name, i) => [String(i), name])];
      field.addEventListener("value-commit", (event) => {
        event.stopPropagation(); const value = event.detail.value;
        if (!this.#setValues([cluster.uid], index, value === "" ? null : this.#domain === null ? Number(value) : this.#domain[Number(value)])) event.preventDefault();
      });
      row.append(cell);
    });
    if (extra) row.append(this.#element("empty-cell"));
    return this.#updateRow(row, cluster);
  }

  #updateRow(row, cluster) {
    const label = cluster.name || cluster.uid, name = row.querySelector("th"); name.textContent = label; name.title = cluster.uid;
    row.querySelector("input").setAttribute("aria-label", `Select ${label}`);
    [...row.querySelectorAll("clinker-editable-field")].forEach((field, index) => {
      const value = this.#panel.values[cluster.uid]?.[index] ?? null;
      field.setAttribute("label", `${label}: ${this.#panel.columns[index]}`);
      field.value = this.#domain === null ? value : value === null ? "" : String(this.#domain.indexOf(value));
    });
    return row;
  }
}

function defineAnnotationTable() {
  defineEditableField(); defineSelectionBar(); defineActionMenu(); defineDataTable(); defineAnnotationValueActions();
  if (!customElements.get("clinker-annotation-table")) customElements.define("clinker-annotation-table", AnnotationTable);
}

/** Shared CSV/TSV/JSON boundary for creation and replacement. */
function parseAnnotationTable(text, categorical = false) {
  const raw = text.trim();
  if (raw.startsWith("{")) {
    const { columns, values } = JSON.parse(raw);
    return { columns, values };
  }
  let parse = csvParseRows;
  for (const token of raw.matchAll(/"(?:[^"]|"")*"|([,\t\r\n])/g)) {
    if (token[1]) { parse = token[1] === "\t" ? tsvParseRows : csvParseRows; break; }
  }
  const rows = parse(raw);
  if (rows.length < 2 || rows[0].length < 2) throw new Error("Include a header and at least one cluster row.");
  const columns = rows.shift().slice(1).map((label) => label.trim()), values = Object.create(null);
  for (const row of rows) {
    if (row.every((cell) => !cell.trim())) continue;
    const uid = row[0].trim();
    if (Object.hasOwn(values, uid)) throw new Error(`Duplicate cluster ID '${uid}'.`);
    if (row.length > columns.length + 1) throw new Error(`Too many values for '${uid}'.`);
    values[uid] = columns.map((_, index) => {
      const value = (row[index + 1] ?? "").trim();
      if (!value) return null;
      if (categorical) return value;
      const numeric = Number(value);
      if (!Number.isFinite(numeric)) throw new Error(`Expected a numeric value for '${uid}'.`);
      return numeric;
    });
  }
  return { columns, values };
}

function createAnnotationPanel({ type, title, scale, data }, panels = []) {
  const raw = data.trim();
  const contents = type === "groups" ? { groups: JSON.parse(raw).groups } : parseAnnotationTable(raw, scale === "categorical");
  const ids = new Set(panels.map(({ uid }) => uid));
  let serial = 1;
  while (ids.has(`${type}-${serial}`)) serial++;
  return { uid: `${type}-${serial}`, type,
    title: title.trim() || (type === "groups" ? "Groups" : type === "bar" ? "Bar chart" : "Heatmap"),
    ...contents, ...(type === "groups" ? {} : { scale: { type: scale } }) };
}

function annotationTableText(panel, clusters) {
  return csvFormatRows([["cluster", ...panel.columns], ...clusters.map((cluster) =>
    [cluster.uid, ...(panel.values[cluster.uid] ?? panel.columns.map(() => null))])]);
}

function mapAnnotationValues(values, map) {
  return Object.fromEntries(Object.entries(values).map(([uid, row]) => [uid, row.map(map)]));
}

const [markup, css] = await Promise.all([
  "<template id=\"editor\">\n<div class=\"cm-annotation-layout\">\n  <clinker-annotation-toolbar>\n    <clinker-search-field label=\"Filter records\" placeholder=\"Find a cluster\"></clinker-search-field>\n  </clinker-annotation-toolbar>\n  <div class=\"cm-annotation-content\">\n    <div class=\"cm-data-error\" data-data-error role=\"alert\"></div>\n    <div class=\"cm-annotation-data-body\"></div>\n  </div>\n</div>\n</template>\n\n<template id=\"matrix\">\n<clinker-annotation-table></clinker-annotation-table>\n</template>\n\n<template id=\"categories\">\n<clinker-category-editor slot=\"tools\" class=\"cm-annotation-categories\"></clinker-category-editor>\n</template>\n\n<template id=\"creation\">\n<clinker-annotation-data-form></clinker-annotation-data-form>\n</template>\n\n<template id=\"replacement\">\n<clinker-annotation-data-form slot=\"tools\" mode=\"replace\" hidden></clinker-annotation-data-form>\n</template>\n\n<template id=\"groups\">\n<div class=\"cm-data-panel cm-annotation-group-list\" role=\"region\" aria-label=\"Editable annotation data\" tabindex=\"-1\"><clinker-annotation-group-editor></clinker-annotation-group-editor></div>\n</template>\n\n<template id=\"empty\">\n<p>No annotation panels. Add a panel to attach data to cluster rows.</p>\n</template>\n",
  ".cm-annotation-data-editor { --_cm-control-height:30px; --_cm-control-padding:4px 8px; --_cm-control-radius:6px; outline:none; display:flex; flex-direction:column; gap:8px; min-height:0; }\n.cm-annotation-layout, .cm-annotation-content { display:flex; flex-direction:column; flex:1; gap:8px; min-width:0; min-height:0; }\nclinker-annotation-editor[embedded], clinker-annotation-editor[embedded] > .cm-annotation-layout { display:contents; }\nclinker-annotation-editor[hidden] { display:none; }\nclinker-annotation-editor[embedded] .cm-annotation-layout > clinker-annotation-toolbar { grid-area:tools; min-width:0; }\nclinker-annotation-editor[embedded] .cm-annotation-content { grid-area:content; }\n.cm-annotation-data-body { display:flex; flex-direction:column; flex:1; min-height:0; overflow:hidden; }\n.cm-annotation-group-list { overflow:auto; }\n.cm-annotation-data-editor [data-data-error] { color:var(--cm-danger,var(--_cm-danger,#b42318)); }\n.cm-annotation-data-editor [data-data-error]:empty { display:none; }\n",
]);
const template = createTemplate(markup);
const HTMLElementBase = globalThis.HTMLElement || class {};

class AnnotationEditor extends HTMLElementBase {
  #state = { data: null, clusters: [], patch: () => {} };
  #templates;
  #toolbar;
  #body;
  #search;
  #status;
  #applying = false;
  #selected = null;
  #creating = false;
  #view;
  #viewKey;
  #categories;
  #form;
  #run = createActionRunner((message) => { this.#status.textContent = message; });

  connectedCallback() {
    const root = this.getRootNode();
    installDataEditorStyle(root);
    this.classList.add("cm-ui", "cm-data-ui", "cm-annotation-data-editor");
    installStyle(root, "clinker-annotation-editor-style", css);
    if (!this.#body) this.#mount();
    this.#render();
  }
  disconnectedCallback() { this.#cancelEdit(); }
  set state(value) { this.#state = value; if (this.#body && !this.#applying) this.#render(); }
  get value() { return this.#selected; }
  select(uid) { this.#selected = uid; this.#creating = false; if (this.#body) this.#render(); }
  clearSelection() { this.#view?.clearSelection?.(); }
  reset() {
    this.#cancelEdit(); this.#selected = this.#viewKey = null;
    this.#state = { data: null, clusters: [], patch: () => {} };
    this.#creating = false;
    if (this.#search) this.#search.value = "";
    if (this.#body) this.#render();
  }
  get #palette() { return this.#state.annotationApi.ANNOTATION_PALETTES.Qualitative; }
  #scale(panel) { return this.#state.annotationApi.resolveAnnotationScale(panel); }
  #panels() { return this.#state.data?.annotationPanels || []; }
  #current() { return this.#panels().find((panel) => panel.uid === this.#selected); }
  #element(name) { return this.ownerDocument.importNode(this.#templates.getElementById(name).content, true).firstElementChild; }
  #cancelEdit() { this.#view?.cancelEdit(); this.#categories?.cancelEdit(); }
  #commit(operations) {
    this.#applying = true;
    try { return this.#run(() => { this.#state.patch(operations); this.#render(); }); }
    finally { this.#applying = false; }
  }
  #update(changes) { return this.#commit([{ type: "annotationPanels.update", ids: [this.#selected], changes }]); }
  #categoryScale(panel, domain, colours) {
    return { ...panel.scale, type: "categorical", domain, colours: colours.length ? colours : [...this.#palette] };
  }
  #mount() {
    this.#templates = template(this.ownerDocument);
    this.append(this.#element("editor"));
    this.#toolbar = this.querySelector("clinker-annotation-toolbar");
    this.#body = this.querySelector(".cm-annotation-data-body");
    this.#status = this.querySelector("[data-data-error]");
    this.#search = this.querySelector("clinker-search-field");
    this.#search.addEventListener("input", () => this.#syncView());
    this.#toolbar.addEventListener("panel-change", (event) => this.select(event.detail.value));
    this.#toolbar.addEventListener("panel-action", (event) => this.#action(event.detail.action));
    this.#body.addEventListener("columns-change", (event) => { if (this.#update(event.detail) === false) event.preventDefault(); });
    this.#body.addEventListener("rows-change", (event) => {
      const panel = this.#current(), categorical = panel.scale?.type === "categorical", scale = this.#scale(panel);
      if (this.#commit([
        ...(categorical && !panel.scale?.domain ? [{ type: "annotationPanels.update", ids: [this.#selected],
          changes: { scale: this.#categoryScale(panel, scale.domain, scale.domain.map(scale.colour)) } }] : []),
        { type: "annotationPanels.updateRows", uid: this.#selected, values: event.detail.values },
      ]) === false) event.preventDefault();
    });
    this.#body.addEventListener("value-commit", (event) => {
      let changes;
      if (event.target === this.#categories) {
        const { value, replacements } = event.detail, mapping = new Map(replacements);
        changes = { scale: this.#categoryScale(this.#current(), value.map(({ name }) => name), value.map(({ colour }) => colour)) };
        if (mapping.size) changes.values = mapAnnotationValues(this.#current().values, (item) => mapping.has(item) ? mapping.get(item) : item);
      } else if (event.target === this.#view && this.#current()?.type === "groups") changes = { groups: event.detail.value };
      else return;
      if (this.#update(changes) === false) event.preventDefault();
    });
    this.#body.addEventListener("data-submit", (event) => { if (event.target === this.#form) this.#submit(event.detail); });
    this.#body.addEventListener("data-cancel", (event) => { if (event.target === this.#form) this.#form.hidden = true; });
  }
  #action(action) {
    if (action === "add") { this.#creating = !this.#creating; this.#render(); }
    else if (action === "replace" && this.#form) {
      this.#form.hidden = !this.#form.hidden;
      if (!this.#form.hidden) {
        this.#form.value = { data: annotationTableText(this.#current(), this.#state.clusters) };
        this.#form.error = ""; this.#form.focus();
      }
    } else if (action === "show") this.#update({ show: this.#current().show === false });
    else if (action === "remove") this.#commit([{ type: "annotationPanels.delete", ids: [this.#selected] }]);
    else if (action === "left" || action === "right") {
      const ids = this.#panels().map((panel) => panel.uid), index = ids.indexOf(this.#selected), next = index + (action === "left" ? -1 : 1);
      if (index < 0 || next < 0 || next >= ids.length) return;
      [ids[index], ids[next]] = [ids[next], ids[index]];
      this.#commit([{ type: "annotationPanels.reorder", ids }]);
    }
  }
  #submit(detail) {
    try {
      if (this.#creating) {
        const panel = createAnnotationPanel(detail, this.#panels());
        this.#state.patch([{ type: "annotationPanels.create", panel }]); this.select(panel.uid); return;
      }
      const source = this.#current(), categorical = source.scale?.type === "categorical", changes = parseAnnotationTable(detail.data, categorical);
      if (categorical) {
        const scale = this.#scale(source), domain = new Set(scale.domain);
        for (const row of Object.values(changes.values)) for (const value of row) if (value !== null) domain.add(value);
        const names = [...domain]; changes.scale = this.#categoryScale(source, names, names.map((value, i) =>
          scale.domain.includes(value) ? scale.colour(value) : this.#palette[i % this.#palette.length]));
      }
      if (this.#update(changes)) this.#form.hidden = true;
    } catch (exception) { this.#form.error = exception.message; }
  }
  #syncView() {
    if (!this.#view) return;
    const panel = this.#current(), clusters = this.#state.clusters;
    this.#view.state = panel.type === "groups" ? { groups: panel.groups, colour: panel.colour, clusters, tree: this.#state.data.tree }
      : { panel, clusters, filter: this.#search.value, domain: panel.scale?.type === "categorical" ? this.#scale(panel).domain : null };
  }
  #render() {
    const panels = this.#panels();
    if (!panels.some((panel) => panel.uid === this.#selected)) this.#selected = panels[0]?.uid ?? null;
    const panel = this.#current();
    this.#toolbar.state = { panels, value: this.#selected, creating: this.#creating };
    this.#search.hidden = !panel || this.#creating || panel.type === "groups";
    const categorical = panel?.scale?.type === "categorical";
    const key = this.#creating ? "create" : panel ? JSON.stringify([panel.uid, panel.type]) : "empty";
    if (key !== this.#viewKey) {
      this.#cancelEdit(); this.#body.replaceChildren(); this.#view = this.#categories = this.#form = null; this.#viewKey = key;
      if (this.#creating) {
        this.#form = this.#element("creation"); this.#form.clusters = this.#state.clusters; this.#body.append(this.#form);
      } else if (!panel) this.#body.append(this.#element("empty"));
      else if (panel.type === "groups") {
        const region = this.#element("groups"); this.#body.append(region); this.#view = region.firstElementChild;
      } else {
        this.#view = this.#element("matrix"); this.#body.append(this.#view);
        this.#form = this.#element("replacement"); this.#view.append(this.#form);
      }
    }
    this.#syncView();
    if (this.#view && panel.type !== "groups" && categorical && !this.#categories) {
      this.#categories = this.#element("categories"); this.#view.prepend(this.#categories);
      this.#categories.palette = this.#palette;
    } else if (!categorical && this.#categories) {
      this.#categories.cancelEdit(); this.#categories.remove(); this.#categories = null;
    }
    if (this.#categories) {
      const { domain, colour } = this.#scale(panel);
      this.#categories.value = domain.map((name) => ({ name, colour: colour(name) }));
    }
  }
}

function defineAnnotationEditor() {
  defineSearchField(); defineCategoryEditor(); defineAnnotationToolbar(); defineAnnotationDataForm(); defineAnnotationGroupEditor(); defineAnnotationTable();
  if (!customElements.get("clinker-annotation-editor")) customElements.define("clinker-annotation-editor", AnnotationEditor);
}

export { DATA_EDITOR_CSS as D, RowSelection as R, installStyle as a, defineDataTable as b, defineEditableField as c, defineActionMenu as d, defineSelectionBar as e, defineAnnotationEditor as f, DATA_ROW_HEIGHT as g, defineAnnotationDataForm as h, installDataEditorStyle as i, createAnnotationPanel as j };
