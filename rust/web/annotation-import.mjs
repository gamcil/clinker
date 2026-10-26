const NUMBER = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;
const SKIP = "__skip__";

export function parseAnnotationTable(text, delimiter) {
  const source = String(text).replace(/^\uFEFF/, "");
  if (delimiter !== "," && delimiter !== "\t") throw new Error("Choose comma or tab as the separator.");
  const rows = [];
  let row = [], field = "", quoted = false, closedQuote = false;
  for (let i = 0; i < source.length; i += 1) {
    const char = source[i];
    if (quoted) {
      if (char === '"' && source[i + 1] === '"') { field += '"'; i += 1; }
      else if (char === '"') { quoted = false; closedQuote = true; }
      else field += char;
      continue;
    }
    if (char === '"') {
      if (field.length || closedQuote) throw new Error(`Unexpected quote near character ${i + 1}.`);
      quoted = true;
    } else if (char === delimiter) {
      row.push(field); field = ""; closedQuote = false;
    } else if (char === "\n" || char === "\r") {
      if (char === "\r" && source[i + 1] === "\n") i += 1;
      row.push(field); field = ""; closedQuote = false;
      if (row.some(value => value.trim() !== "")) rows.push(row);
      row = [];
    } else {
      if (closedQuote && !/\s/.test(char)) throw new Error(`Unexpected text after a quoted field near character ${i + 1}.`);
      if (!closedQuote) field += char;
    }
  }
  if (quoted) throw new Error("The table contains an unclosed quoted field.");
  row.push(field);
  if (row.some(value => value.trim() !== "")) rows.push(row);
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
  if (text !== undefined) node.textContent = text;
  return node;
}

// Preview and import use the same plan, including scale inference on matched rows.
export function planAnnotationPanels({ table, columns, assignments, mode = "tracks", matrixScale = "auto", title = "" }) {
  if (!columns.length) return [];
  const groups = mode === "matrix" ? [columns] : columns.map(column => [column]);
  return groups.map(group => {
    const rawValues = [...assignments.keys()].flatMap(index => group.map(column => table.records[index][column.index])).filter(value => value !== "");
    const numeric = rawValues.length > 0 && rawValues.every(value => NUMBER.test(value) && Number.isFinite(Number(value)));
    const binary = numeric && rawValues.every(value => Number(value) === 0 || Number(value) === 1);
    let scaleKind = mode === "matrix" ? matrixScale : group[0].valueKind;
    if (scaleKind === "auto") {
      const kinds = new Set(group.map(column => column.valueKind));
      if (kinds.size > 1) throw new Error("Select numeric columns for a measurement matrix, or choose Categories to treat all selected values as labels.");
      scaleKind = binary ? "binary" : kinds.has("continuous") ? "continuous" : "categorical";
    }
    if (scaleKind === "continuous" && rawValues.length && !numeric) throw new Error("A numeric heatmap can only include numeric values. Skip text columns or choose Categories.");
    if (scaleKind === "binary" && rawValues.length && !binary) throw new Error("Presence / absence requires values of 0 or 1. Blank cells remain missing.");
    const scale = scaleKind === "binary"
      ? { type: "categorical", domain: ["0", "1"], colours: ["#f7fbff", "#08306b"] }
      : { type: scaleKind };
    const values = Object.fromEntries([...assignments].map(([index, uid]) => [uid, group.map(column => {
      const raw = table.records[index][column.index];
      return raw === "" ? null : scale.type === "continuous" ? Number(raw) : scaleKind === "binary" ? String(Number(raw)) : raw;
    })]));
    const type = mode === "matrix" ? "heatmap" : group[0].kind;
    const present = Object.values(values).flat().filter(value => value !== null);
    let key, scaleSummary;
    if (scale.type === "categorical") {
      const categories = scale.domain || [...new Set(present)];
      key = scaleKind === "binary" ? "0 = absent; 1 = present" : categories.length ? `${categories.slice(0, 5).join(", ")}${categories.length > 5 ? ` … (${categories.length} categories)` : ""}` : "No matched values";
      scaleSummary = scaleKind === "binary" ? "0 / 1" : `${categories.length} categor${categories.length === 1 ? "y" : "ies"}`;
    } else if (present.length) {
      let min = type === "bar" ? 0 : Infinity, max = type === "bar" ? 0 : -Infinity;
      for (const value of present) { min = Math.min(min, value); max = Math.max(max, value); }
      if (min === max) {
        if (type === "bar") { min = 0; max = 1; }
        else { min -= 0.5; max += 0.5; }
      }
      key = `${min} to ${max}`;
    } else key = "No non-missing values";
    return {
      // A track is one variable; a matrix has a title plus per-variable labels.
      panel: { type, title: mode === "matrix" ? title.trim() : group[0].name,
        showColumnLabels: mode === "matrix", columns: group.map(column => column.name), values, scale },
      label: mode === "matrix" ? title.trim() || "Heatmap matrix" : group[0].name,
      scaleLabel: scaleKind === "binary" ? "Presence / absence" : scale.type === "continuous" ? "Numeric" : "Categories",
      key,
      scaleSummary: scaleSummary || key,
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
  const modeInput = dialog.querySelector("#annotation-mode");
  const matrixOptions = dialog.querySelector("#annotation-matrix-options");
  const matrixScale = dialog.querySelector("#annotation-matrix-scale");
  const matrixTitle = dialog.querySelector("#annotation-matrix-title");
  const columnNote = dialog.querySelector("#annotation-column-note");
  const preview = dialog.querySelector("#annotation-preview");
  const previewSummary = dialog.querySelector("#annotation-preview-summary");
  const parsingOptions = dialog.querySelector("#annotation-parsing-options");
  const apply = dialog.querySelector("#apply-annotations");
  let table = null, clusters = [], matchEntries = [], file = null;
  let readVersion = 0;
  const columnChoices = { tracks: new Map(), matrix: new Map() };

  const selectedColumns = () => [...columnsBody.querySelectorAll("tr")].map(row => ({
    name: table.headers[Number(row.dataset.index)], index: Number(row.dataset.index),
    kind: row.querySelector("select").value, valueKind: row.querySelector("select").dataset.kind,
  })).filter(column => column.kind !== "skip");
  const makePlan = assignments => planAnnotationPanels({ table, columns: selectedColumns(), assignments,
    mode: modeInput.value, matrixScale: matrixScale.value, title: matrixTitle.value });

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
    const scaleCells = new Map([...columnsBody.querySelectorAll("tr")].map(row => [table?.headers[Number(row.dataset.index)], row.querySelector(".annotation-column-scale")]));
    for (const cell of scaleCells.values()) { cell.textContent = "—"; cell.removeAttribute("title"); }
    let plan = [], invalidPlan = false;
    if (table) {
      try {
        plan = makePlan(assignments);
        previewSummary.textContent = plan.length
          ? modeInput.value === "matrix" ? `${plan[0].panel.columns.length} columns · Shared ${plan[0].scaleLabel.toLowerCase()} scale: ${plan[0].scaleSummary}. No normalisation.`
            : `${plan.length} track${plan.length === 1 ? "" : "s"} with independent scales.`
          : "Select at least one column to import.";
        for (const item of plan) {
          for (const name of item.panel.columns) {
            const cell = scaleCells.get(name);
            cell.textContent = modeInput.value === "matrix" ? "Shared" : item.scaleSummary;
            cell.title = `${item.scaleLabel}: ${item.key}`;
          }
        }
        if (unresolved || duplicateTargets.size) previewSummary.textContent += " Resolve row matching to import.";
      } catch (error) { invalidPlan = true; previewSummary.textContent = error.message; }
    }
    previewSummary.dataset.tone = invalidPlan ? "error" : "info";
    apply.textContent = modeInput.value === "matrix" ? "Import heatmap" : plan.length ? `Import ${plan.length} track${plan.length === 1 ? "" : "s"}` : "Import tracks";
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

  const renderColumns = () => {
    columnsBody.replaceChildren();
    const keyIndex = Number(keyInput.value);
    for (let index = 0; index < table.headers.length; index += 1) {
      if (index === keyIndex) continue;
      const kind = columnKind(table.records, index);
      const row = element(document, "tr");
      row.dataset.index = String(index);
      row.append(element(document, "td", "annotation-column-name", table.headers[index]));
      const cell = element(document, "td");
      const select = element(document, "select", "annotation-column-kind");
      select.dataset.kind = kind;
      select.setAttribute("aria-label", `Import ${table.headers[index]} as`);
      const choices = modeInput.value === "matrix" ? [["heatmap", "Include"], ["skip", "Skip"]]
        : [...(kind === "continuous" ? [["bar", "Bar chart"]] : []), ["heatmap", "Heatmap"], ["skip", "Skip"]];
      for (const [value, label] of choices) {
        const option = element(document, "option", "", label); option.value = value; select.append(option);
      }
      if (skipColumnByDefault(table.headers[index], table.records, index, kind) || modeInput.value === "matrix" && kind !== "continuous") select.value = "skip";
      if (columnChoices[modeInput.value].has(index)) select.value = columnChoices[modeInput.value].get(index);
      select.addEventListener("change", () => {
        columnChoices[modeInput.value].set(index, select.value);
        updateSummary();
      });
      cell.append(select); row.append(cell, element(document, "td", "annotation-column-scale", "—")); columnsBody.append(row);
    }
    columnsWrap.hidden = false;
    updateSummary();
  };

  const loadFile = async () => {
    const version = ++readVersion;
    file = fileInput.files?.[0] || null;
    table = null; matchEntries = [];
    columnsBody.replaceChildren();
    columnsWrap.hidden = true; reviewWrap.hidden = true; preview.hidden = true;
    keyInput.replaceChildren(); keyInput.disabled = true;
    columnChoices.tracks.clear(); columnChoices.matrix.clear();
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
      matchEntries = matchRows(table.records, Number(keyInput.value), clusters);
      parsingOptions.open = matchEntries.some(entry => !entry.cluster);
      renderMatches(); renderColumns();
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
  modeInput.addEventListener("change", () => {
    matrixOptions.hidden = modeInput.value !== "matrix";
    columnNote.textContent = modeInput.value === "matrix"
      ? "Use columns with comparable units. One shared colour scale."
      : "Each column has its own scale.";
    if (table) renderColumns();
  });
  matrixScale.addEventListener("change", updateSummary);
  matrixTitle.addEventListener("input", updateSummary);
  keyInput.addEventListener("change", () => {
    if (!table) return;
    matchEntries = matchRows(table.records, Number(keyInput.value), clusters);
    renderMatches(); renderColumns();
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
    modeInput.value = "tracks"; matrixOptions.hidden = true;
    matrixOptions.open = false; parsingOptions.open = false;
    matrixScale.value = "auto"; matrixTitle.value = "";
    columnChoices.tracks.clear(); columnChoices.matrix.clear();
    columnNote.textContent = "Each column has its own scale.";
    reviewBody.replaceChildren();
    reviewWrap.hidden = true;
    summary.textContent = "Choose a TSV or CSV file to begin.";
    table = null; matchEntries = []; file = null;
    apply.disabled = true;
    apply.textContent = "Import tracks";
    button.closest('details')?.querySelector('summary').focus();
  });
  return { open: () => button.click(), close: () => dialog.close() };
}
