import { parseDelimitedRows } from "./delimited.mjs?v=bf5a7880d0c1d0f2";

const NUMBER = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
const SKIP = "__skip__";

export function parseAnnotationTable(text, delimiter) {
  const rows = parseDelimitedRows(text, delimiter);
  if (rows.length < 2) throw new Error("The table needs a header row and at least one data row.");
  const headers = rows[0].map(value => value.trim());
  if (headers.some(value => !value)) throw new Error("Every column needs a non-empty header.");
  const normalized = headers.map(value => value.toLocaleLowerCase());
  if (new Set(normalized).size !== headers.length) throw new Error("Column headers must be unique (ignoring capitalization).");
  const records = rows.slice(1).map((values, index) => {
    if (values.length !== headers.length) throw new Error(`Row ${index + 2} has ${values.length} fields; expected ${headers.length}.`);
    return values.map(value => value.trim());
  });
  return { headers, records };
}

function columnKind(records, index) {
  const values = records.map(row => row[index]).filter(value => value !== "");
  return values.length && values.every(value => NUMBER.test(value) && Number.isFinite(Number(value)))
    ? "continuous" : "categorical";
}

function skipColumnByDefault(header, records, index, kind) {
  if (/(^|[_\s-])(id|name|url|uri|file|accession)([_\s-]|$)/i.test(header)) return true;
  if (kind !== "categorical") return false;
  const values = records.map(row => row[index]).filter(value => value !== "");
  const unique = new Set(values).size;
  return unique > 20 && unique / values.length > 0.5;
}

function aliasesFor(cluster) {
  const aliases = [cluster.uid, cluster.name, cluster.source_file];
  if (cluster.source_file) {
    const basename = cluster.source_file.split(/[\\/]/).at(-1);
    aliases.push(basename, basename.replace(/\.[^.]+$/, ""));
  }
  for (const locus of cluster.loci || []) aliases.push(locus.name);
  return [...new Set(aliases.filter(value => typeof value === "string" && value.trim()))];
}

function matchRows(records, keyIndex, clusters) {
  const aliases = new Map();
  for (const cluster of clusters) for (const alias of aliasesFor(cluster)) {
    const key = alias.trim().toLocaleLowerCase();
    if (!aliases.has(key)) aliases.set(key, new Set());
    aliases.get(key).add(cluster.uid);
  }
  const entries = records.map((values, index) => {
    const key = values[keyIndex];
    const candidates = key ? [...(aliases.get(key.toLocaleLowerCase()) || [])] : [];
    return { index, key, candidates, cluster: candidates.length === 1 ? candidates[0] : "" };
  });
  const byCluster = new Map();
  for (const entry of entries) if (entry.cluster) {
    if (!byCluster.has(entry.cluster)) byCluster.set(entry.cluster, []);
    byCluster.get(entry.cluster).push(entry);
  }
  for (const duplicates of byCluster.values()) if (duplicates.length > 1) {
    for (const entry of duplicates) { entry.cluster = ""; entry.candidates = [duplicates[0].cluster || entry.candidates[0]]; entry.conflict = true; }
  }
  return entries;
}

function element(document, tag, className, text) {
  const node = document.createElement(tag);
  if (className) node.className = className;
  if (["input", "select", "textarea"].includes(tag)) node.classList.add("app-field");
  if (tag === "select") node.classList.add("flat-control", "compact-control");
  if (tag === "button") node.classList.add("app-button");
  if (text !== undefined) node.textContent = text;
  return node;
}

const HEATMAP_PALETTES = {
  blue: ["#f7fbff", "#08306b"],
  diverging: ["#2166ac", "#f7f7f7", "#b2182b"],
  viridis: ["#440154", "#21918c", "#fde725"],
  qualitative: ["#4477aa", "#ee6677", "#228833", "#ccbb44", "#66ccee", "#aa3377", "#bbbbbb"],
  muted: ["#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677", "#882255"],
};

// Preview and import use the same panel plan and matched-row scale inference.
export function planAnnotationPanels({ table, columns, assignments, heatmapGroups = [] }) {
  const groupById = new Map(heatmapGroups.map(group => [group.id, group]));
  const outputs = [];
  const grouped = new Map();
  for (const column of columns) {
    if (column.display === "skip") continue;
    const groupId = column.display.startsWith("group:") ? column.display.slice(6) : null;
    if (groupId) {
      const config = groupById.get(groupId);
      if (!config) continue;
      if (!grouped.has(groupId)) {
        const entry = { columns: [], config };
        grouped.set(groupId, entry); outputs.push(entry);
      }
      grouped.get(groupId).columns.push(column);
    } else outputs.push({ columns: [column], config: null, type: column.display });
  }

  return outputs.map(({ columns: group, config, type: separateType }) => {
    const heatmap = Boolean(config) || separateType === "heatmap";
    const bar = separateType === "bar";
    const rawValues = [...assignments.keys()].flatMap(index => group.map(column => table.records[index][column.index])).filter(value => value !== "");
    const numeric = rawValues.length > 0 && rawValues.every(value => NUMBER.test(value) && Number.isFinite(Number(value)));
    const binary = numeric && rawValues.every(value => Number(value) === 0 || Number(value) === 1);
    let scaleKind = config?.scale || group[0].valueKind;
    if (scaleKind === "auto") {
      const kinds = new Set(group.map(column => column.valueKind));
      if (kinds.size > 1) throw new Error(`Heatmap group “${config.title}” mixes numeric and text columns. Use one value type per group.`);
      scaleKind = binary ? "binary" : kinds.has("continuous") ? "continuous" : "categorical";
    }
    if (scaleKind === "continuous" && rawValues.length && !numeric) throw new Error(`Heatmap group “${config?.title || group[0].name}” contains non-numeric values. Skip text columns or choose Categories.`);
    if (scaleKind === "binary" && rawValues.length && !binary) throw new Error("Presence / absence requires values of 0 or 1. Blank cells remain missing.");
    const scale = scaleKind === "binary"
      ? { type: "categorical", domain: ["0", "1"], colours: HEATMAP_PALETTES.blue }
      : { type: scaleKind };
    const paletteName = config?.palette ?? group[0].palette;
    if (paletteName && paletteName !== "default") scale.colours = HEATMAP_PALETTES[paletteName];
    const values = Object.fromEntries([...assignments].map(([index, uid]) => [uid, group.map(column => {
      const raw = table.records[index][column.index];
      return raw === "" ? null : scale.type === "continuous" ? Number(raw) : scaleKind === "binary" ? String(Number(raw)) : raw;
    })]));
    const type = heatmap ? "heatmap" : "bar";
    const present = Object.values(values).flat().filter(value => value !== null);
    if (paletteName === "diverging" && scale.type === "continuous" && present.length) {
      const bound = Math.max(0.5, ...present.map(value => Math.abs(value)));
      scale.domain = [-bound, bound];
    }
    let key, scaleSummary;
    if (scale.type === "categorical") {
      const categories = scale.domain || [...new Set(present)];
      key = scaleKind === "binary" ? "0 = absent; 1 = present" : categories.length ? `${categories.slice(0, 5).join(", ")}${categories.length > 5 ? ` … (${categories.length} categories)` : ""}` : "No matched values";
      scaleSummary = scaleKind === "binary" ? "0 / 1" : `${categories.length} categor${categories.length === 1 ? "y" : "ies"}`;
    } else if (present.length) {
      let min = bar ? 0 : Infinity, max = bar ? 0 : -Infinity;
      for (const value of present) { min = Math.min(min, value); max = Math.max(max, value); }
      if (min === max) {
        if (bar) { min = 0; max = 1; }
        else { min -= 0.5; max += 0.5; }
      }
      key = `${min} to ${max}`;
    } else key = "No non-missing values";
    if (scale.type === "continuous" && scale.domain) key = `${scale.domain[0]} to ${scale.domain[1]}`;
    const title = config ? config.title.trim() : group[0].name;
    const panel = { type, title, showColumnLabels: Boolean(config), columns: group.map(column => column.name), values, scale };
    if (!config && group[0].colour) panel.colour = group[0].colour;
    return {
      panel, label: title, grouped: Boolean(config), columnIndices: group.map(column => column.index),
      scaleLabel: scaleKind === "binary" ? "Presence / absence" : scale.type === "continuous" ? "Numeric" : "Categories",
      key, scaleSummary: scaleSummary || key,
    };
  });
}

export function mountAnnotationImport({ button, getChart, onImported = () => {} }) {
  const dialog = document.querySelector("#annotation-dialog");
  const fileInput = dialog.querySelector("#annotation-file");
  const chooseFileButton = dialog.querySelector("#choose-annotation-file");
  const fileName = dialog.querySelector("#annotation-file-name");
  const delimiterInput = dialog.querySelector("#annotation-delimiter");
  const keyInput = dialog.querySelector("#annotation-key-column");
  const summary = dialog.querySelector("#annotation-match-summary");
  const reviewBody = dialog.querySelector("#annotation-match-rows");
  const reviewWrap = dialog.querySelector("#annotation-match-review");
  const columnsBody = dialog.querySelector("#annotation-columns");
  const columnsWrap = dialog.querySelector("#annotation-column-settings");
  const groupSelectedButton = dialog.querySelector("#group-selected-columns");
  const separateSelectedButton = dialog.querySelector("#separate-selected-columns");
  const groupTargetInput = dialog.querySelector("#annotation-group-target");
  const columnActions = dialog.querySelector("#annotation-column-actions");
  const selectionCount = dialog.querySelector("#annotation-selection-count");
  const clearSelectionButton = dialog.querySelector("#clear-annotation-selection");
  const separatePanels = dialog.querySelector("#annotation-separate-panels");
  const preview = dialog.querySelector("#annotation-preview");
  const previewSummary = dialog.querySelector("#annotation-preview-summary");
  const parsingOptions = dialog.querySelector("#annotation-parsing-options");
  const apply = dialog.querySelector("#apply-annotations");
  let table = null, clusters = [], matchEntries = [], file = null;
  let readVersion = 0;
  let currentKeyIndex = null;
  const columnDisplays = new Map();
  const columnGroups = new Map();
  const columnNames = new Map();
  const columnColours = new Map();
  const columnPalettes = new Map();
  const selectedColumnIndices = new Set();
  const heatmapGroups = [];
  const panelOrder = [];
  const collapsedGroups = new Set();
  let selectionAnchor = null;
  let nextGroupId = 1;
  const panelColours = HEATMAP_PALETTES.muted;
  const defaultPanelColour = index => panelColours[index % panelColours.length];

  const columnDefinition = index => {
    if (index === Number(keyInput.value)) return null;
    const valueKind = columnKind(table.records, index);
    const groupId = columnGroups.get(index);
    const display = groupId ? `group:${groupId}` : columnDisplays.get(index) ?? defaultDisplay(index, valueKind);
    const palette = columnPalettes.get(index) ?? (valueKind === "continuous" ? "viridis" : "qualitative");
    return { name: (columnNames.get(index) ?? table.headers[index]).trim(), colour: columnColours.get(index) ?? defaultPanelColour(index), palette, index, display, valueKind };
  };
  const selectedColumns = () => panelOrder.flatMap(key => {
    if (key.startsWith("column:")) {
      const column = columnDefinition(Number(key.slice("column:".length)));
      return column ? [column] : [];
    }
    const id = key.slice("group:".length);
    return table.headers.flatMap((_, index) => columnGroups.get(index) === id ? [columnDefinition(index)] : []);
  });
  const makePlan = assignments => planAnnotationPanels({ table, columns: selectedColumns(), assignments, heatmapGroups });
  const columnOrderKey = index => `column:${index}`;
  const groupOrderKey = id => `group:${id}`;
  const groupMembers = id => table.headers.flatMap((_, index) => columnGroups.get(index) === id && index !== Number(keyInput.value) ? [index] : []);
  const syncPanelOrder = () => {
    if (!table) { panelOrder.length = 0; return; }
    const valid = new Set();
    for (const group of heatmapGroups) if (groupMembers(group.id).length) valid.add(groupOrderKey(group.id));
    table.headers.forEach((_, index) => {
      if (index !== Number(keyInput.value) && !columnGroups.has(index)) valid.add(columnOrderKey(index));
    });
    const seen = new Set();
    for (let index = panelOrder.length - 1; index >= 0; index -= 1) {
      if (!valid.has(panelOrder[index]) || seen.has(panelOrder[index])) panelOrder.splice(index, 1);
      else seen.add(panelOrder[index]);
    }
    for (const group of heatmapGroups) {
      const key = groupOrderKey(group.id);
      if (valid.has(key) && !seen.has(key)) { panelOrder.push(key); seen.add(key); }
    }
    table.headers.forEach((_, index) => {
      const key = columnOrderKey(index);
      if (valid.has(key) && !seen.has(key)) { panelOrder.push(key); seen.add(key); }
    });
  };

  const addHeatmapGroup = ({ title, scale = "auto", palette = "default" } = {}) => {
    const group = { id: `heatmap-${nextGroupId++}`, title: title || `Heatmap ${nextGroupId - 1}`, scale, palette };
    heatmapGroups.push(group);
    return group;
  };

  const selectedAssignments = () => new Map(matchEntries.filter(entry => entry.cluster).map(entry => [entry.index, entry.cluster]));
  const clusterOptions = (select, selected) => {
    select.append(element(document, "option", "", "Choose a cluster…"));
    const skip = element(document, "option", "", "Skip this row"); skip.value = SKIP; select.append(skip);
    for (const cluster of clusters) {
      const label = `${cluster.name || cluster.uid} · ${cluster.uid}${cluster.source_file ? ` · ${cluster.source_file}` : ""}`;
      const option = element(document, "option", "", label); option.value = cluster.uid; select.append(option);
    }
    select.value = selected || "";
  };

  const updateSummary = () => {
    const assignments = selectedAssignments();
    const targets = [...assignments.values()];
    const duplicateTargets = new Set(targets.filter((uid, index) => targets.indexOf(uid) !== index));
    const unresolved = matchEntries.filter(entry => !entry.cluster && entry.choice !== SKIP).length;
    const skipped = matchEntries.filter(entry => entry.choice === SKIP).length;
    const matched = assignments.size - duplicateTargets.size;
    const missing = clusters.length - matched;
    summary.textContent = `${matched} of ${matchEntries.length} rows matched.`;
    if (missing) summary.textContent += ` ${missing} clusters missing.`;
    if (skipped) summary.textContent += ` ${skipped} rows skipped.`;
    if (unresolved || duplicateTargets.size) {
      summary.textContent += unresolved ? ` Resolve or skip ${unresolved} unmatched/ambiguous row${unresolved === 1 ? "" : "s"}.` : " A cluster is assigned more than once; change or skip duplicate rows.";
    }
    preview.hidden = !table;
    const scaleCells = new Map([...columnsWrap.querySelectorAll("[data-index]")].map(row => [Number(row.dataset.index), row.querySelector(".annotation-column-scale")]));
    for (const cell of scaleCells.values()) { cell.textContent = "—"; cell.removeAttribute("title"); }
    let plan = [], invalidPlan = false;
    if (table) {
      try {
        const activeColumns = selectedColumns().filter(column => column.display !== "skip");
        const names = activeColumns.map(column => column.name);
        if (names.some(name => !name)) throw new Error("Column names cannot be empty.");
        if (new Set(names.map(name => name.toLocaleLowerCase())).size !== names.length) throw new Error("Imported column names must be unique.");
        const usedGroups = new Set(activeColumns.map(column => column.display.startsWith("group:") ? column.display.slice(6) : null).filter(Boolean));
        const groupTitles = heatmapGroups.filter(group => usedGroups.has(group.id)).map(group => group.title.trim());
        if (groupTitles.some(title => !title)) throw new Error("Heatmap group names cannot be empty.");
        if (new Set(groupTitles.map(title => title.toLocaleLowerCase())).size !== groupTitles.length) throw new Error("Heatmap group names must be unique.");
        plan = makePlan(assignments);
        const groupedPanels = plan.filter(item => item.grouped);
        const separatePanels = plan.length - groupedPanels.length;
        const groupedColumns = groupedPanels.reduce((total, item) => total + item.panel.columns.length, 0);
        previewSummary.textContent = plan.length
          ? `${plan.length} panel${plan.length === 1 ? "" : "s"}: ${separatePanels} separate, ${groupedPanels.length} heatmap group${groupedPanels.length === 1 ? "" : "s"}${groupedColumns ? ` (${groupedColumns} columns)` : ""}.`
          : "Choose a display for at least one column.";
        for (const item of plan) {
          for (const index of item.columnIndices) {
            const cell = scaleCells.get(index);
            if (!cell) continue;
            cell.textContent = item.grouped ? "Shared" : item.scaleSummary;
            cell.title = `${item.scaleLabel}: ${item.key}`;
          }
        }
        if (unresolved || duplicateTargets.size) previewSummary.textContent += " Resolve row matching to import.";
      } catch (error) { invalidPlan = true; previewSummary.textContent = error.message; }
    }
    previewSummary.dataset.tone = invalidPlan ? "error" : "info";
    apply.textContent = plan.length ? `Import ${plan.length} panel${plan.length === 1 ? "" : "s"}` : "Import panels";
    apply.disabled = !table || !plan.length || !assignments.size || invalidPlan || unresolved > 0 || duplicateTargets.size > 0;
    return { assignments, unresolved, duplicateTargets };
  };

  const renderMatches = () => {
    reviewBody.replaceChildren();
    const needsReview = matchEntries.filter(entry => entry.candidates.length !== 1 || entry.conflict);
    reviewWrap.hidden = needsReview.length === 0;
    for (const entry of needsReview) {
      const row = element(document, "tr");
      row.append(element(document, "td", "annotation-match-key", entry.key || "(blank)"));
      const cell = element(document, "td");
      const select = element(document, "select", "annotation-match-select");
      select.setAttribute("aria-label", `Match table key ${entry.key || "blank"} to a plotted cluster`);
      const defaultChoice = entry.choice || entry.cluster || "";
      clusterOptions(select, defaultChoice);
      entry.choice = select.value;
      select.addEventListener("change", () => { entry.choice = select.value; entry.cluster = select.value === SKIP ? "" : select.value; updateSummary(); });
      cell.append(select); row.append(cell); reviewBody.append(row);
    }
    updateSummary();
  };

  const defaultDisplay = (index, kind) => {
    if (skipColumnByDefault(table.headers[index], table.records, index, kind)) return "skip";
    return kind === "continuous" ? "bar" : "heatmap";
  };

  const paletteOptions = kind => kind === "continuous"
    ? [["default", "Default"], ["viridis", "Viridis"], ["diverging", "Blue–red (zero-centred)"]]
    : [["default", "Default"], ["qualitative", "Qualitative"], ["muted", "Muted"]];

  const updateSelectionActions = () => {
    const selected = [...selectedColumnIndices];
    const selectedGroups = [...new Set(selected.map(index => columnGroups.get(index)).filter(Boolean))];
    const hasGroupedSelection = selectedGroups.length > 0;
    const showActions = selected.length >= 2 || hasGroupedSelection;
    columnActions.hidden = !showActions;
    selectionCount.textContent = `${selected.length} ${selected.length === 1 ? "column" : "columns"} selected`;
    const previousTarget = groupTargetInput.value;
    const suggestedTarget = selectedGroups.length === 1 ? selectedGroups[0] : "";
    groupTargetInput.replaceChildren();
    const newGroup = element(document, "option", "", heatmapGroups.length ? "New heatmap group" : "Create heatmap group");
    newGroup.value = ""; groupTargetInput.append(newGroup);
    for (const group of heatmapGroups) {
      const option = element(document, "option", "", group.title);
      option.value = group.id; groupTargetInput.append(option);
    }
    groupTargetInput.value = heatmapGroups.some(group => group.id === suggestedTarget)
      ? suggestedTarget
      : heatmapGroups.some(group => group.id === previousTarget) ? previousTarget : "";
    groupTargetInput.hidden = heatmapGroups.length === 0 || selected.length < 2;
    const targetGroup = groupTargetInput.value;
    const selectionNeedsGrouping = selected.length >= 2 && (!targetGroup || selected.some(index => columnGroups.get(index) !== targetGroup));
    groupSelectedButton.textContent = targetGroup ? "Add to heatmap" : "Group as heatmap";
    groupSelectedButton.hidden = !selectionNeedsGrouping;
    separateSelectedButton.hidden = !hasGroupedSelection;
    for (const checkbox of columnsBody.querySelectorAll("[data-group-select]")) {
      const members = table.headers.flatMap((_, index) => columnGroups.get(index) === checkbox.dataset.groupSelect ? [index] : []);
      const selectedMembers = members.filter(index => selectedColumnIndices.has(index)).length;
      checkbox.checked = members.length > 0 && selectedMembers === members.length;
      checkbox.indeterminate = selectedMembers > 0 && selectedMembers < members.length;
      checkbox.closest("tr")?.classList.toggle("is-selected", members.length > 0 && selectedMembers === members.length);
    }
  };

  const createColumnNameField = index => {
    const field = element(document, "clinker-editable-field", "annotation-column-name-input");
    field.setAttribute("compact", ""); field.setAttribute("required", "");
    field.setAttribute("label", `Display name for ${table.headers[index]}`);
    field.value = columnNames.get(index) ?? table.headers[index];
    field.addEventListener("value-commit", event => {
      const name = event.detail.value.trim();
      columnNames.set(index, name); field.value = name; updateSummary();
    });
    return field;
  };

  const createPanelOrderHandle = (key, label) => {
    const handle = element(document, "button", "app-button quiet-action compact-control annotation-panel-reorder", "⠿");
    handle.type = "button"; handle.draggable = true; handle.dataset.reorderPanel = ""; handle.dataset.orderKey = key;
    handle.setAttribute("aria-label", `Reorder ${label}; use Arrow Up or Arrow Down to move it`);
    handle.title = "Drag to reorder; Arrow Up or Arrow Down moves one place";
    return handle;
  };

  const createColumnCheckbox = index => {
    const checkbox = element(document, "input", "annotation-column-select");
    checkbox.type = "checkbox";
    checkbox.checked = selectedColumnIndices.has(index);
    checkbox.setAttribute("aria-label", `Select ${table.headers[index]}`);
    checkbox.addEventListener("click", event => {
      const orderedIndices = [...columnsBody.querySelectorAll("tr[data-index]")].map(row => Number(row.dataset.index));
      const anchorPosition = orderedIndices.indexOf(selectionAnchor);
      const targetPosition = orderedIndices.indexOf(index);
      if (event.shiftKey && anchorPosition >= 0 && targetPosition >= 0) {
        for (const selectedIndex of orderedIndices.slice(Math.min(anchorPosition, targetPosition), Math.max(anchorPosition, targetPosition) + 1)) {
          if (checkbox.checked) selectedColumnIndices.add(selectedIndex);
          else selectedColumnIndices.delete(selectedIndex);
        }
      } else {
        if (checkbox.checked) selectedColumnIndices.add(index);
        else selectedColumnIndices.delete(index);
        selectionAnchor = index;
      }
      for (const row of columnsBody.querySelectorAll("tr[data-index]")) {
        const rowIndex = Number(row.dataset.index);
        row.classList.toggle("is-selected", selectedColumnIndices.has(rowIndex));
        const rowCheckbox = row.querySelector(".annotation-column-select");
        if (rowCheckbox) rowCheckbox.checked = selectedColumnIndices.has(rowIndex);
      }
      updateSelectionActions();
    });
    return checkbox;
  };

  const renderColumns = () => {
    columnsBody.replaceChildren();
    const appendGroup = group => {
      const indices = table.headers.flatMap((_, index) => columnGroups.get(index) === group.id && index !== Number(keyInput.value) ? [index] : []);
      if (!indices.length) return;
      const row = element(document, "tr", "annotation-panel-group-row"); row.dataset.groupId = group.id;
      row.dataset.orderKey = groupOrderKey(group.id); row.dataset.panelRow = "";
      const selectCell = element(document, "td", "annotation-column-select-cell");
      const groupSelect = element(document, "input", "annotation-column-select");
      groupSelect.type = "checkbox"; groupSelect.dataset.groupSelect = group.id;
      groupSelect.setAttribute("aria-label", `Select all columns in ${group.title}`);
      groupSelect.addEventListener("change", () => {
        for (const index of indices) {
          if (groupSelect.checked) selectedColumnIndices.add(index);
          else selectedColumnIndices.delete(index);
        }
        selectionAnchor = indices.at(-1) ?? selectionAnchor;
        renderColumns();
      });
      const reorderHandle = createPanelOrderHandle(groupOrderKey(group.id), group.title);
      selectCell.append(groupSelect, reorderHandle);
      const toggleGroup = element(document, "button", "app-button quiet-action compact-control annotation-group-toggle");
      toggleGroup.type = "button";
      toggleGroup.title = collapsedGroups.has(group.id) ? "Expand heatmap panel" : "Collapse heatmap panel";
      toggleGroup.setAttribute("aria-label", toggleGroup.title);
      toggleGroup.setAttribute("aria-expanded", String(!collapsedGroups.has(group.id)));
      const disclosureIcon = document.createElementNS("http://www.w3.org/2000/svg", "svg");
      disclosureIcon.setAttribute("viewBox", "0 0 16 16"); disclosureIcon.setAttribute("aria-hidden", "true");
      const disclosurePath = document.createElementNS("http://www.w3.org/2000/svg", "path");
      disclosurePath.setAttribute("d", "m5 3 5 5-5 5"); disclosureIcon.append(disclosurePath); toggleGroup.append(disclosureIcon);
      toggleGroup.addEventListener("click", () => {
        if (collapsedGroups.has(group.id)) collapsedGroups.delete(group.id);
        else collapsedGroups.add(group.id);
        renderColumns();
      });
      selectCell.append(toggleGroup);

      const nameCell = element(document, "td", "annotation-panel-name");
      const name = element(document, "div", "annotation-panel-name-content");
      const title = element(document, "clinker-editable-field", "annotation-heatmap-group-title");
      title.setAttribute("compact", ""); title.setAttribute("required", "");
      title.setAttribute("label", "Heatmap panel name"); title.value = group.title;
      title.addEventListener("value-commit", event => {
        group.title = event.detail.value.trim(); title.value = group.title;
        reorderHandle.setAttribute("aria-label", `Reorder ${group.title}; use Arrow Up or Arrow Down to move it`);
        updateSelectionActions(); updateSummary();
      });
      name.append(title); nameCell.append(name);

      const displayCell = element(document, "td", "annotation-panel-type", "Heatmap");
      const paletteCell = element(document, "td", "annotation-column-colour-cell");
      const palette = element(document, "select", "annotation-group-palette");
      palette.setAttribute("aria-label", `Palette for ${group.title}`);
      const groupKinds = new Set(indices.map(index => columnKind(table.records, index)));
      const paletteKind = group.scale === "continuous" ? "continuous" : group.scale === "categorical" || group.scale === "binary"
        ? "categorical" : groupKinds.size === 1 ? [...groupKinds][0] : "categorical";
      for (const [value, label] of paletteOptions(paletteKind)) {
        const option = element(document, "option", "", label); option.value = value; palette.append(option);
      }
      palette.value = group.palette;
      palette.addEventListener("change", () => { group.palette = palette.value; updateSummary(); });
      paletteCell.append(palette);
      const scaleCell = element(document, "td", "annotation-column-scale");
      const scale = element(document, "select", "annotation-group-scale");
      scale.setAttribute("aria-label", `Scale for ${group.title}`);
      for (const [value, label] of [["auto", "Automatic"], ["continuous", "Numeric"], ["binary", "Presence / absence"], ["categorical", "Categories"]]) {
        const option = element(document, "option", "", label); option.value = value; scale.append(option);
      }
      scale.value = group.scale;
      scale.addEventListener("change", () => {
        group.scale = scale.value;
        const nextKind = group.scale === "continuous" ? "continuous" : group.scale === "categorical" || group.scale === "binary"
          ? "categorical" : paletteKind;
        if (!paletteOptions(nextKind).some(([value]) => value === group.palette)) group.palette = nextKind === "continuous" ? "viridis" : "qualitative";
        renderColumns();
      });
      const scaleControls = element(document, "div", "annotation-group-scale-controls");
      const removeGroup = element(document, "button", "app-button icon-button annotation-remove-group", "×");
      removeGroup.type = "button"; removeGroup.title = "Remove heatmap panel"; removeGroup.setAttribute("aria-label", "Remove heatmap panel");
      removeGroup.addEventListener("click", () => {
        const orderIndex = panelOrder.indexOf(groupOrderKey(group.id));
        if (orderIndex >= 0) panelOrder.splice(orderIndex, 1, ...indices.map(columnOrderKey));
        for (const index of indices) {
          columnGroups.delete(index);
          columnDisplays.set(index, defaultDisplay(index, columnKind(table.records, index)));
        }
        heatmapGroups.splice(heatmapGroups.indexOf(group), 1);
        collapsedGroups.delete(group.id);
        selectionAnchor = null;
        renderAllColumns();
      });
      scaleControls.append(scale, removeGroup); scaleCell.append(scaleControls);
      row.append(selectCell, nameCell, displayCell, paletteCell, scaleCell);
      row.classList.toggle("is-selected", indices.length > 0 && indices.every(index => selectedColumnIndices.has(index)));
      columnsBody.append(row);
      if (!collapsedGroups.has(group.id)) for (const index of indices) appendColumn(index, group.id);
    };

    const appendColumn = (index, groupId = null) => {
      const kind = columnKind(table.records, index);
      const row = element(document, "tr", groupId ? "annotation-heatmap-member" : ""); row.dataset.index = String(index);
      row.dataset.orderKey = groupId ? groupOrderKey(groupId) : columnOrderKey(index);
      if (!groupId) row.dataset.panelRow = "";
      row.classList.toggle("is-selected", selectedColumnIndices.has(index));
      const checkCell = element(document, "td", "annotation-column-select-cell");
      checkCell.append(createColumnCheckbox(index));
      if (!groupId) checkCell.append(createPanelOrderHandle(columnOrderKey(index), table.headers[index]));
      const nameCell = element(document, "td", groupId ? "annotation-column-name annotation-group-member-name" : "annotation-column-name");
      if (groupId) nameCell.append(createColumnNameField(index));
      else {
        const nameContent = element(document, "div", "annotation-panel-name-content");
        nameContent.append(createColumnNameField(index));
        nameCell.append(nameContent);
      }
      const displayCell = element(document, "td");
      const colourCell = element(document, "td", "annotation-column-colour-cell");
      const scaleCell = element(document, "td", "annotation-column-scale", groupId ? "Shared" : "—");
      if (groupId) {
        scaleCell.title = "Shared heatmap scale";
      } else {
        const select = element(document, "select", "annotation-column-kind"); select.dataset.kind = kind;
        select.setAttribute("aria-label", `Display ${table.headers[index]} as`);
        for (const [value, label] of [["skip", "Skip"], ...(kind === "continuous" ? [["bar", "Bar chart"]] : []), ["heatmap", "Heatmap"]]) {
          const option = element(document, "option", "", label); option.value = value; select.append(option);
        }
        const display = columnDisplays.get(index) ?? defaultDisplay(index, kind);
        columnDisplays.set(index, display); select.value = display;
        const renderColourControl = () => {
          colourCell.replaceChildren();
          if (select.value === "skip") {
            colourCell.append(element(document, "span", "annotation-column-scale", "—"));
          } else if (select.value === "bar") {
            const colourInput = element(document, "input", "annotation-column-colour");
            colourInput.type = "color"; colourInput.value = columnColours.get(index) || defaultPanelColour(index);
            colourInput.setAttribute("aria-label", `Panel colour for ${table.headers[index]}`);
            colourInput.addEventListener("input", () => { columnColours.set(index, colourInput.value); updateSummary(); });
            colourCell.append(colourInput);
          } else {
            const palette = element(document, "select", "annotation-group-palette");
            palette.setAttribute("aria-label", `Palette for ${table.headers[index]}`);
            for (const [value, label] of paletteOptions(kind)) {
              const option = element(document, "option", "", label); option.value = value; palette.append(option);
            }
            palette.value = columnPalettes.get(index) || (kind === "continuous" ? "viridis" : "qualitative");
            palette.addEventListener("change", () => { columnPalettes.set(index, palette.value); updateSummary(); });
            colourCell.append(palette);
          }
        };
        select.addEventListener("change", () => {
          columnDisplays.set(index, select.value); renderColourControl(); updateSummary();
        });
        displayCell.append(select); renderColourControl();
      }
      row.append(checkCell, nameCell, displayCell, colourCell, scaleCell); columnsBody.append(row);
    };

    syncPanelOrder();
    for (const key of panelOrder) {
      if (key.startsWith("group:")) {
        const group = heatmapGroups.find(item => groupOrderKey(item.id) === key);
        if (group) appendGroup(group);
      } else {
        const index = Number(key.slice("column:".length));
        if (index !== Number(keyInput.value) && !columnGroups.has(index)) appendColumn(index);
      }
    }
    const hasColumns = table.headers.some((_, index) => index !== Number(keyInput.value));
    separatePanels.hidden = !hasColumns && heatmapGroups.length === 0;
    columnsWrap.hidden = false;
    updateSelectionActions(); updateSummary();
  };

  let draggedPanelKey = null;
  let dropPanelKey = null;
  let dropAfter = false;
  const clearDropTarget = () => {
    for (const row of columnsBody.querySelectorAll(".is-drop-before, .is-drop-after")) row.classList.remove("is-drop-before", "is-drop-after");
    dropPanelKey = null;
  };
  const movePanel = (sourceKey, targetKey, after, restoreFocus = false) => {
    if (sourceKey === targetKey) return;
    const sourceIndex = panelOrder.indexOf(sourceKey);
    if (sourceIndex < 0 || !panelOrder.includes(targetKey)) return;
    panelOrder.splice(sourceIndex, 1);
    const targetIndex = panelOrder.indexOf(targetKey);
    panelOrder.splice(targetIndex + (after ? 1 : 0), 0, sourceKey);
    renderColumns();
    if (restoreFocus) {
      const handle = [...columnsBody.querySelectorAll("[data-reorder-panel]")].find(item => item.dataset.orderKey === sourceKey);
      handle?.focus({ preventScroll: true });
    }
  };
  const movePanelBy = (key, direction) => {
    const index = panelOrder.indexOf(key), target = index + direction;
    if (index < 0 || target < 0 || target >= panelOrder.length) return;
    panelOrder.splice(index, 1); panelOrder.splice(target, 0, key);
    renderColumns();
    const handle = [...columnsBody.querySelectorAll("[data-reorder-panel]")].find(item => item.dataset.orderKey === key);
    handle?.focus({ preventScroll: true });
  };
  const rowsForOrderKey = key => [...columnsBody.querySelectorAll("tr[data-order-key]")].filter(row => row.dataset.orderKey === key);
  const panelRowFor = key => [...columnsBody.querySelectorAll("tr[data-panel-row]")].find(row => row.dataset.orderKey === key);
  const updateDropTarget = (key, clientY) => {
    clearDropTarget();
    const targetRow = panelRowFor(key);
    if (!targetRow) return;
    const rows = rowsForOrderKey(key), first = rows[0]?.getBoundingClientRect(), last = rows.at(-1)?.getBoundingClientRect();
    const midpoint = first && last ? (first.top + last.bottom) / 2 : (targetRow.getBoundingClientRect().top + targetRow.getBoundingClientRect().bottom) / 2;
    dropAfter = clientY >= midpoint; dropPanelKey = key;
    targetRow.classList.add(dropAfter ? "is-drop-after" : "is-drop-before");
  };

  columnsBody.addEventListener("keydown", event => {
    const handle = event.target.closest("[data-reorder-panel]");
    if (!handle || !["ArrowUp", "ArrowDown"].includes(event.key)) return;
    event.preventDefault(); movePanelBy(handle.dataset.orderKey, event.key === "ArrowUp" ? -1 : 1);
  });
  columnsBody.addEventListener("dragstart", event => {
    const handle = event.target.closest("[data-reorder-panel]");
    if (!handle) return;
    draggedPanelKey = handle.dataset.orderKey;
    event.dataTransfer.effectAllowed = "move";
    event.dataTransfer.setData("text/plain", draggedPanelKey);
  });
  columnsBody.addEventListener("dragover", event => {
    if (!draggedPanelKey) return;
    const row = event.target.closest("tr[data-order-key]");
    if (!row) return;
    const key = row.dataset.orderKey;
    if (key === draggedPanelKey) { clearDropTarget(); return; }
    event.preventDefault(); event.dataTransfer.dropEffect = "move"; updateDropTarget(key, event.clientY);
  });
  columnsBody.addEventListener("drop", event => {
    if (draggedPanelKey && dropPanelKey) {
      event.preventDefault(); movePanel(draggedPanelKey, dropPanelKey, dropAfter);
    }
    draggedPanelKey = null; clearDropTarget();
  });
  columnsBody.addEventListener("dragend", () => { draggedPanelKey = null; clearDropTarget(); });

  const renderAllColumns = () => renderColumns();

  const pruneEmptyGroups = () => {
    for (let index = heatmapGroups.length - 1; index >= 0; index -= 1) {
      const group = heatmapGroups[index];
      if (![...columnGroups.values()].includes(group.id)) {
        heatmapGroups.splice(index, 1);
      }
    }
  };

  const moveSelectedToGroup = (groupId, placeAtSelection = false) => {
    const selected = new Set(selectedColumnIndices);
    let precedingKeys = [];
    if (placeAtSelection) {
      const sourceIndex = panelOrder.findIndex(key => key.startsWith("column:")
        ? selected.has(Number(key.slice("column:".length)))
        : groupMembers(key.slice("group:".length)).some(index => selected.has(index)));
      let boundary = sourceIndex < 0 ? panelOrder.length : sourceIndex;
      const sourceKey = panelOrder[sourceIndex];
      if (sourceKey?.startsWith("group:")) {
        const sourceMembers = groupMembers(sourceKey.slice("group:".length));
        if (sourceMembers.some(index => !selected.has(index))) boundary += 1;
      }
      precedingKeys = panelOrder.slice(0, boundary).filter(key => {
        if (key.startsWith("column:")) return !selected.has(Number(key.slice("column:".length)));
        const id = key.slice("group:".length);
        return id === groupId || groupMembers(id).some(index => !selected.has(index));
      });
    }
    for (const index of selected) {
      if (index === Number(keyInput.value)) continue;
      columnGroups.set(index, groupId); columnDisplays.delete(index);
    }
    selectedColumnIndices.clear();
    selectionAnchor = null;
    pruneEmptyGroups(); syncPanelOrder();
    if (placeAtSelection) {
      const key = groupOrderKey(groupId);
      const currentIndex = panelOrder.indexOf(key);
      if (currentIndex >= 0) panelOrder.splice(currentIndex, 1);
      const insertionIndex = precedingKeys.filter(item => panelOrder.includes(item)).length;
      panelOrder.splice(insertionIndex, 0, key);
    }
    renderAllColumns();
  };

  const groupSelectedColumns = () => {
    if (!selectedColumnIndices.size) return;
    if (heatmapGroups.some(group => group.id === groupTargetInput.value)) {
      moveSelectedToGroup(groupTargetInput.value);
      return;
    }
    const group = addHeatmapGroup();
    moveSelectedToGroup(group.id, true);
  };

  const separateSelectedColumns = () => {
    const selected = new Set([...selectedColumnIndices].filter(index => columnGroups.has(index)));
    const nextOrder = [];
    for (const key of panelOrder) {
      if (!key.startsWith("group:")) { nextOrder.push(key); continue; }
      const id = key.slice("group:".length), members = groupMembers(id);
      const separating = members.filter(index => selected.has(index));
      if (!separating.length) { nextOrder.push(key); continue; }
      const remaining = members.filter(index => !selected.has(index));
      if (remaining.length) nextOrder.push(key);
      nextOrder.push(...separating.map(columnOrderKey));
    }
    panelOrder.splice(0, panelOrder.length, ...nextOrder);
    for (const index of selected) {
      columnGroups.delete(index);
      columnDisplays.set(index, defaultDisplay(index, columnKind(table.records, index), "tracks"));
    }
    selectedColumnIndices.clear(); selectionAnchor = null; pruneEmptyGroups(); syncPanelOrder(); renderAllColumns();
  };

  const loadFile = async () => {
    const version = ++readVersion;
    file = fileInput.files?.[0] || null;
    table = null; matchEntries = []; currentKeyIndex = null;
    columnsBody.replaceChildren();
    columnsWrap.hidden = true; reviewWrap.hidden = true; preview.hidden = true;
    keyInput.replaceChildren(); keyInput.disabled = true;
    columnDisplays.clear(); columnGroups.clear(); selectedColumnIndices.clear();
    columnNames.clear(); columnColours.clear(); columnPalettes.clear(); panelOrder.length = 0;
    heatmapGroups.length = 0; collapsedGroups.clear(); selectionAnchor = null; nextGroupId = 1;
    updateSelectionActions();
    apply.disabled = true;
    if (!file) { fileName.textContent = "None selected"; keyInput.disabled = true; summary.textContent = "Choose a TSV or CSV file to begin."; return; }
    fileName.textContent = file.name;
    try {
      const delimiter = delimiterInput.value === "auto" ? (/\.tsv$/i.test(file.name) ? "\t" : ",") : delimiterInput.value;
      const text = await file.text();
      if (version !== readVersion) return;
      table = parseAnnotationTable(text, delimiter);
      const chart = getChart();
      if (!chart?.data()) throw new Error("Run an analysis or open a project before importing annotations.");
      clusters = chart.data().clusters;
      keyInput.replaceChildren();
      table.headers.forEach((header, index) => {
        const option = element(document, "option", "", header); option.value = String(index); keyInput.append(option);
      });
      keyInput.disabled = false;
      const preferred = table.headers.findIndex(header => /^(cluster(_?id)?|cluster name|file(name)?|source_file|record(_?id)?|id)$/i.test(header));
      keyInput.value = String(preferred < 0 ? 0 : preferred);
      currentKeyIndex = Number(keyInput.value);
      table.headers.forEach((_, index) => { if (index !== currentKeyIndex) panelOrder.push(columnOrderKey(index)); });
      matchEntries = matchRows(table.records, Number(keyInput.value), clusters);
      parsingOptions.open = matchEntries.some(entry => !entry.cluster);
      renderMatches(); renderAllColumns();
    } catch (error) {
      if (version !== readVersion) return;
      table = null;
      columnsWrap.hidden = true; reviewWrap.hidden = true; preview.hidden = true;
      apply.disabled = true; keyInput.disabled = true;
      summary.textContent = error.message || String(error);
    }
  };

  button.addEventListener("click", () => {
    if (!getChart()?.data()) return;
    clusters = getChart().data().clusters;
    dialog.showModal();
  });
  fileInput.addEventListener("change", loadFile);
  chooseFileButton.addEventListener("click", () => fileInput.click());
  delimiterInput.addEventListener("change", loadFile);
  groupSelectedButton.addEventListener("click", groupSelectedColumns);
  clearSelectionButton.addEventListener("click", () => {
    selectedColumnIndices.clear(); selectionAnchor = null; renderAllColumns();
  });
  groupTargetInput.addEventListener("change", () => {
    const target = groupTargetInput.value;
    groupSelectedButton.textContent = target ? "Add to heatmap" : "Group as heatmap";
    groupSelectedButton.hidden = selectedColumnIndices.size < 2 || (target && [...selectedColumnIndices].every(index => columnGroups.get(index) === target));
  });
  separateSelectedButton.addEventListener("click", separateSelectedColumns);
  keyInput.addEventListener("change", () => {
    if (!table) return;
    const nextKeyIndex = Number(keyInput.value);
    if (currentKeyIndex !== null && nextKeyIndex !== currentKeyIndex) {
      if (!columnGroups.has(currentKeyIndex) && !columnDisplays.has(currentKeyIndex)) {
        const display = defaultDisplay(currentKeyIndex, columnKind(table.records, currentKeyIndex));
        if (display.startsWith("group:")) columnGroups.set(currentKeyIndex, display.slice(6));
        else columnDisplays.set(currentKeyIndex, display);
      }
      columnGroups.delete(nextKeyIndex); columnDisplays.delete(nextKeyIndex); selectedColumnIndices.delete(nextKeyIndex);
      currentKeyIndex = nextKeyIndex; pruneEmptyGroups();
      syncPanelOrder();
    }
    matchEntries = matchRows(table.records, Number(keyInput.value), clusters);
    renderMatches(); renderAllColumns();
  });
  apply.addEventListener("click", () => {
    try {
      const { assignments, unresolved, duplicateTargets } = updateSummary();
      if (!table || unresolved || duplicateTargets.size) return;
      if (apply.disabled) return;
      const plan = makePlan(assignments);
      const existing = new Set((getChart().data().annotationPanels || []).map(panel => panel.uid));
      const nextUid = index => {
        // Keep IDs readable, including untitled matrices in the editor selector.
        const base = plan[index].label;
        let uid = base, suffix = 2;
        while (existing.has(uid)) uid = `${base} (${suffix++})`;
        existing.add(uid);
        return uid;
      };
      const panels = plan.map(({ panel }, index) => ({ ...panel, uid: nextUid(index) }));
      const chart = getChart();
      chart.patch(panels.map(panel => ({ type: "annotationPanels.create", panel })));
      chart.config({ annotations: { enabled: true } });
      const importedFileName = file.name;
      dialog.close();
      onImported(`${panels.length} annotation panel${panels.length === 1 ? "" : "s"} imported from ${importedFileName}.`);
    } catch (error) {
      summary.textContent = `Could not import annotations: ${error.message || String(error)}`;
    }
  });

  dialog.addEventListener("close", () => {
    readVersion += 1;
    fileInput.value = "";
    fileName.textContent = "None selected";
    keyInput.replaceChildren();
    keyInput.disabled = true;
    columnsBody.replaceChildren();
    columnsWrap.hidden = true;
    preview.hidden = true;
    parsingOptions.open = false;
    columnDisplays.clear(); columnGroups.clear(); selectedColumnIndices.clear();
    columnNames.clear(); columnColours.clear(); columnPalettes.clear(); panelOrder.length = 0;
    heatmapGroups.length = 0; collapsedGroups.clear(); selectionAnchor = null; nextGroupId = 1;
    updateSelectionActions();
    reviewBody.replaceChildren();
    reviewWrap.hidden = true;
    summary.textContent = "Choose a TSV or CSV file to begin.";
    table = null; matchEntries = []; file = null; currentKeyIndex = null;
    apply.disabled = true;
    apply.textContent = "Import panels";
    button.closest('details')?.querySelector('summary').focus();
  });
  return { open: () => button.click(), close: () => dialog.close() };
}
