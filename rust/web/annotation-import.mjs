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
  const apply = dialog.querySelector("#apply-annotations");
  let table = null, clusters = [], matchEntries = [], file = null;

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
    summary.textContent = `${matched} of ${matchEntries.length} table rows matched; ${missing} plotted clusters will have missing values; ${skipped} table rows skipped.`;
    if (unresolved || duplicateTargets.size) {
      summary.textContent += unresolved ? ` Resolve or skip ${unresolved} unmatched/ambiguous row${unresolved === 1 ? "" : "s"}.` : " A cluster is assigned more than once; change or skip duplicate rows.";
    }
    const hasSelectedColumn = [...columnsBody.querySelectorAll("select[data-kind]")].some(select => select.value !== "skip");
    apply.disabled = !table || !hasSelectedColumn || unresolved > 0 || duplicateTargets.size > 0;
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
      const defaultChoice = entry.choice || (entry.candidates.length === 1 ? entry.candidates[0] : "");
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
      row.append(element(document, "td", "annotation-column-name", table.headers[index]));
      row.append(element(document, "td", "", kind === "continuous" ? "Numeric" : "Categorical"));
      const cell = element(document, "td");
      const select = element(document, "select", "annotation-column-kind");
      select.dataset.kind = kind;
      select.setAttribute("aria-label", `Import ${table.headers[index]} as`);
      for (const [value, label] of [["heatmap", "Heatmap"], ...(kind === "continuous" ? [["bar", "Bar chart"]] : []), ["skip", "Skip"]]) {
        const option = element(document, "option", "", label); option.value = value; select.append(option);
      }
      if (skipColumnByDefault(table.headers[index], table.records, index, kind)) select.value = "skip";
      select.addEventListener("change", updateSummary);
      cell.append(select); row.append(cell); columnsBody.append(row);
    }
    columnsWrap.hidden = false;
    updateSummary();
  };

  const loadFile = async () => {
    file = fileInput.files?.[0] || null;
    table = null; matchEntries = [];
    columnsWrap.hidden = true; reviewWrap.hidden = true;
    apply.disabled = true;
    if (!file) { fileName.textContent = "None selected"; keyInput.disabled = true; summary.textContent = "Choose a TSV or CSV file to begin."; return; }
    fileName.textContent = file.name;
    try {
      const delimiter = delimiterInput.value === "auto" ? (/\.tsv$/i.test(file.name) ? "\t" : ",") : delimiterInput.value;
      table = parseAnnotationTable(await file.text(), delimiter);
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
      renderMatches(); renderColumns();
    } catch (error) {
      table = null;
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
  keyInput.addEventListener("change", () => {
    if (!table) return;
    matchEntries = matchRows(table.records, Number(keyInput.value), clusters);
    renderMatches();
  });
  apply.addEventListener("click", () => {
    try {
      const { assignments, unresolved, duplicateTargets } = updateSummary();
      if (!table || unresolved || duplicateTargets.size) return;
      const columnKinds = [...columnsBody.querySelectorAll("tr")].map(row => ({
        name: row.querySelector(".annotation-column-name").textContent,
        index: table.headers.indexOf(row.querySelector(".annotation-column-name").textContent),
        kind: row.querySelector("select").value,
        valueKind: row.querySelector("select").dataset.kind,
      })).filter(column => column.kind !== "skip");
      if (!columnKinds.length) return;
      const existing = new Set((getChart().data().annotationPanels || []).map(panel => panel.uid));
      const prefix = `import-${Date.now().toString(36)}`;
      const panels = [];
      const nextUid = index => {
        let uid = `${prefix}-${index}`;
        while (existing.has(uid)) uid += "x";
        existing.add(uid);
        return uid;
      };
      const cellValue = (rowIndex, column) => {
        const raw = table.records[rowIndex][column.index];
        return raw === "" ? null : column.valueKind === "continuous" ? Number(raw) : raw;
      };
      for (const valueKind of ["continuous", "categorical"]) {
        const group = columnKinds.filter(column => column.kind === "heatmap" && column.valueKind === valueKind);
        if (!group.length) continue;
        const values = {};
        for (const [rowIndex, clusterUid] of assignments) {
          values[clusterUid] = group.map(column => cellValue(rowIndex, column));
        }
        panels.push({
          uid: nextUid(panels.length), type: "heatmap",
          title: valueKind === "continuous" ? "Numeric annotations" : "Categorical annotations",
          columns: group.map(column => column.name), values,
          scale: { type: valueKind },
        });
      }
      for (const column of columnKinds.filter(column => column.kind === "bar")) {
        const values = {};
        for (const [rowIndex, clusterUid] of assignments) values[clusterUid] = [cellValue(rowIndex, column)];
        panels.push({
          uid: nextUid(panels.length), type: "bar", title: column.name,
          columns: [column.name], values, scale: { type: "continuous" }, showColumnLabels: false,
        });
      }
      const chart = getChart();
      chart.patch(panels.map(panel => ({ type: "annotationPanels.create", panel })));
      chart.config({ annotations: { enabled: true } });
      dialog.close();
      onImported(`${panels.length} annotation panel${panels.length === 1 ? "" : "s"} imported from ${file.name}.`);
    } catch (error) {
      summary.textContent = `Could not import annotations: ${error.message || String(error)}`;
    }
  });

  dialog.addEventListener("close", () => {
    fileInput.value = "";
    fileName.textContent = "None selected";
    keyInput.replaceChildren();
    keyInput.disabled = true;
    columnsBody.replaceChildren();
    columnsWrap.hidden = true;
    reviewBody.replaceChildren();
    reviewWrap.hidden = true;
    summary.textContent = "Choose a TSV or CSV file to begin.";
    table = null; matchEntries = []; file = null;
    apply.disabled = true;
    button.closest('details')?.querySelector('summary').focus();
  });
  return { open: () => button.click(), close: () => dialog.close() };
}
