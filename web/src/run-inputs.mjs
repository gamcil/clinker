import { inspectInputFiles, isSupportedInputFile } from "./input-files.mjs";
import { inputTreeClusters } from "./viewer/tree.mjs";
import { parseTwoColumnCsv } from "./delimited.mjs";
import { Workflow } from "./workflow.mjs";

export function mountRunInputs(container, { onChange, onInputsChanged, onError }) {
  const document = container.ownerDocument;
  const filesInput = container.querySelector("#files");
  const fileDropZone = container.querySelector("#file-drop-zone");
  const fileSummary = container.querySelector("#file-summary");
  const fileInputStatus = container.querySelector("#file-input-status");
  const fileInputStatusSummary = container.querySelector("#file-input-status-summary");
  const fileInputStatusList = container.querySelector("#file-input-status-list");
  const clearInputFilesButton = container.querySelector("#clear-input-files");
  const identityInput = container.querySelector("#identity");
  const prefilterInput = container.querySelector("#prefilter");
  const noAlignInput = container.querySelector("#no-align");
  const fileOrderInput = container.querySelector("#file-order");
  const kmerSizeInput = container.querySelector("#kmer-size");
  const minSharedKmersInput = container.querySelector("#min-shared-kmers");
  const separateRecordsInput = container.querySelector("#separate-records");
  const keepOriginInput = container.querySelector("#keep-origin");
  const recordRangesElement = container.querySelector("#record-ranges");
  const geneFunctionsInput = container.querySelector("#gene-functions");
  const colourMapInput = container.querySelector("#colour-map");
  const geneFunctionsSummary = container.querySelector("#gene-functions-summary");
  const colourMapSummary = container.querySelector("#colour-map-summary");
  let disabled = false;
  let inputsReady = false;
  let inputInspection = null;
  let inputInspectionPromise = Promise.resolve({ ready: false, rows: [] });
  const inputMetadata = new WeakMap();
  const selectedInputFiles = new Map();

  function syncPrefilterControls() {
    prefilterInput.disabled = noAlignInput.disabled || noAlignInput.checked;
    fileOrderInput.disabled = noAlignInput.disabled || noAlignInput.checked;
    const prefilterDisabled = prefilterInput.disabled || !prefilterInput.checked;
    kmerSizeInput.disabled = prefilterDisabled;
    minSharedKmersInput.disabled = prefilterDisabled;
  }

  async function readGroupDefinitions(workflow) {
    const functionFile = geneFunctionsInput.files[0];
    const colourFile = colourMapInput.files[0];
    if (!functionFile) {
      if (colourFile) throw new Error("Select a gene functions CSV before adding a colour map.");
      return [];
    }
    const functions = parseTwoColumnCsv(await workflow.readFile(functionFile, true), "Gene functions CSV");
    const colours = new Map(colourFile
      ? parseTwoColumnCsv(await workflow.readFile(colourFile, true), "Function colours CSV")
      : []);
    const grouped = new Map();
    for (const [gene, label] of functions) {
      if (!grouped.has(label)) grouped.set(label, []);
      grouped.get(label).push(gene);
    }
    return [...grouped].map(([label, genes]) => ({
      label,
      genes,
      colour: colours.get(label) || null,
    }));
  }

  function syncFileSummary() {
    const count = filesInput.files.length;
    fileSummary.textContent = count === 0
      ? ""
      : `${count} input ${count === 1 ? "file" : "files"} selected`;
    fileSummary.hidden = count === 0;
  }

  function setSelectedInputFiles(files, { replace = false } = {}) {
    if (replace) selectedInputFiles.clear();
    for (const file of files) selectedInputFiles.set(file.name.toLocaleLowerCase(), file);
    const transfer = new DataTransfer();
    for (const file of selectedInputFiles.values()) transfer.items.add(file);
    filesInput.files = transfer.files;
  }

  function renderInputValidation(validation) {
    fileInputStatusList.replaceChildren();
    fileInputStatus.hidden = validation.rows.length === 0;
    if (validation.rows.length === 0) return;
    const errorCount = validation.rows.filter(row => !row.ready).length;
    const fileCount = filesInput.files.length;
    const countLabel = `${fileCount} input ${fileCount === 1 ? "file" : "files"}`;
    fileInputStatusSummary.textContent = errorCount === 0
      ? `${countLabel} ready`
      : `${errorCount} ${errorCount === 1 ? "problem" : "problems"} · ${countLabel}`;
    const orderedRows = [
      ...validation.rows.filter(row => !row.ready),
      ...validation.rows.filter(row => row.ready),
    ];
    for (const row of orderedRows) {
      const item = document.createElement("li");
      item.dataset.state = row.ready ? "ready" : "error";
      const name = document.createElement("span");
      name.className = "file-input-name";
      name.textContent = row.name;
      name.title = row.name;
      const state = document.createElement("span");
      state.className = "file-input-state";
      state.textContent = row.detail;
      const remove = document.createElement("button");
      remove.type = "button";
      remove.className = "app-button remove-input-file quiet-action compact-control";
      remove.textContent = "Remove";
      remove.setAttribute("aria-label", `Remove ${row.name}`);
      remove.disabled = disabled;
      remove.addEventListener("click", () => {
        for (const name of row.files || [row.name]) selectedInputFiles.delete(name.toLocaleLowerCase());
        setSelectedInputFiles([], { replace: false });
        handleFileSelection();
      });
      const description = document.createElement("span");
      description.className = "file-input-description";
      description.textContent = row.description || "";
      const label = document.createElement("div");
      label.className = "file-input-label";
      label.append(name, description);
      item.append(label, state, remove);
      if (row.reason) {
        const reason = document.createElement("span");
        reason.className = "file-input-reason";
        reason.textContent = row.reason;
        item.append(reason);
      }
      fileInputStatusList.append(item);
    }
  }

  function emptyRecordRanges(message) {
    recordRangesElement.replaceChildren();
    const empty = document.createElement("p");
    empty.className = "empty-setting";
    empty.textContent = message;
    recordRangesElement.append(empty);
  }

  function recordCoordinate(kind, value, minimum, maximum, description) {
    const wrapper = document.createElement("label");
    wrapper.className = "record-coordinate";
    const input = document.createElement("input");
    input.className = "app-field compact-control";
    input.type = "number";
    input.min = String(minimum);
    input.max = String(maximum);
    input.step = "1";
    input.value = String(value);
    input.dataset.rangeCoordinate = kind;
    input.setAttribute("aria-label", `${description}: ${kind === "start" ? "Start" : "End"} (bp)`);
    wrapper.append(input);
    return wrapper;
  }

  function renderRecordRanges(filesWithRecords) {
    if (filesWithRecords.length === 0) {
      emptyRecordRanges("Choose input files to trim their sequences.");
      return;
    }
    recordRangesElement.replaceChildren();
    const headings = document.createElement("table");
    headings.className = "app-table sequence-regions sequence-region-headings";
    const head = headings.createTHead().insertRow();
    for (const text of ["File / sequence", "Start (bp)", "End (bp)"]) {
      const cell = document.createElement("th");
      cell.scope = "col";
      cell.textContent = text;
      head.append(cell);
    }
    const header = document.createElement("div");
    header.className = "sequence-region-header";
    header.append(headings);
    const rows = document.createElement("div");
    rows.className = "sequence-region-rows";
    recordRangesElement.append(header, rows);
    for (const { file, records } of filesWithRecords) {
      const table = document.createElement("table");
      table.className = "app-table sequence-regions";
      table.setAttribute("aria-label", file.name);
      const columns = document.createElement("colgroup");
      for (let index = 0; index < 3; index += 1) {
        const column = document.createElement("col");
        if (index > 0) column.className = "sequence-coordinate-column";
        columns.append(column);
      }
      table.append(columns);
      rows.append(table);
      const body = table.createTBody();
      if (records.length > 1) table.classList.add("sequence-file-group");
      if (records.length === 0) {
        const cell = body.insertRow().insertCell();
        cell.colSpan = 3;
        cell.className = "empty-setting";
        cell.textContent = `${file.name}: No sequences found; the complete file will be analysed.`;
        continue;
      }
      for (const record of records) {
        const availableStart = record.start ?? 1;
        const availableEnd = record.end ?? record.length;
        const row = body.insertRow();
        row.dataset.recordRange = "";
        row.dataset.file = file.name;
        row.dataset.record = record.name;
        row.dataset.length = String(record.length);
        row.dataset.availableStart = String(availableStart);
        row.dataset.availableEnd = String(availableEnd);
        const description = `${file.name}: ${record.name}`;
        const identity = document.createElement("th");
        identity.scope = "row";
        identity.title = `${description} (${availableStart.toLocaleString()}–${availableEnd.toLocaleString()} bp)`;
        const fileName = document.createElement("span");
        fileName.className = "sequence-file-name";
        fileName.textContent = file.name;
        if (record !== records[0]) {
          fileName.classList.add("sequence-file-name-placeholder");
          fileName.setAttribute("aria-hidden", "true");
        }
        identity.append(fileName);
        const name = document.createElement("span");
        name.className = "record-range-name";
        name.textContent = record.name;
        identity.append(name);
        row.append(identity);
        row.insertCell().append(recordCoordinate("start", availableStart, availableStart, availableEnd, description));
        row.insertCell().append(recordCoordinate("end", availableEnd, availableStart, availableEnd, description));
      }
    }
  }

  function collectRecordRanges() {
    const selections = new Map();
    for (const row of recordRangesElement.querySelectorAll("[data-record-range]")) {
      const availableStart = Number(row.dataset.availableStart);
      const availableEnd = Number(row.dataset.availableEnd);
      const start = Number(row.querySelector('[data-range-coordinate="start"]').value);
      const end = Number(row.querySelector('[data-range-coordinate="end"]').value);
      const description = `${row.dataset.file}: ${row.dataset.record}`;
      if (!Number.isSafeInteger(start) || !Number.isSafeInteger(end)
          || start < availableStart || end < start || end > availableEnd) {
        throw new Error(`${description} requires whole-number coordinates between ${availableStart.toLocaleString()} and ${availableEnd.toLocaleString()}.`);
      }
      const selection = { record: row.dataset.record, start: start - 1, end };
      if (!selections.has(selection.record)) selections.set(selection.record, []);
      selections.get(selection.record).push({
        range: selection,
        isFullRecord: selection.start === availableStart - 1 && selection.end === availableEnd,
      });
    }
    const ranges = [];
    for (const [record, instances] of selections) {
      if (instances.every(instance => instance.isFullRecord)) continue;
      if (instances.length > 1 && instances.some(instance => instance.isFullRecord)) {
        throw new Error(`Record ${record} occurs in multiple files; adjust every instance to the same range.`);
      }
      const [first, ...rest] = instances;
      if (rest.some(instance => instance.range.start !== first.range.start
          || instance.range.end !== first.range.end)) {
        throw new Error(`Record ${record} occurs in multiple files with different ranges.`);
      }
      ranges.push(first.range);
    }
    return ranges;
  }

  function syncAnnotationFiles() {
    const functionFile = geneFunctionsInput.files[0];
    const colourFile = colourMapInput.files[0];
    geneFunctionsSummary.textContent = functionFile?.name || "None selected";
    colourMapInput.disabled = geneFunctionsInput.disabled || !functionFile;
    colourMapSummary.textContent = colourFile?.name || (functionFile ? "None selected" : "Add functions first");
  }

  function handleFileSelection() {
    syncFileSummary();
    inputInspection?.dispose();
    const inspection = new Workflow();
    inputInspection = inspection;
    const files = [...filesInput.files];
    inputsReady = false;
    onChange();
    fileInputStatusList.replaceChildren();
    fileInputStatus.hidden = files.length === 0;
    if (files.length > 0) fileInputStatusSummary.textContent = "Checking selected inputs…";
    recordRangesElement.setAttribute("aria-busy", "true");
    emptyRecordRanges(files.length ? "Reading annotation records…" : "Choose input files to trim their sequences.");
    inputInspectionPromise = inspectInputFiles(files, {
      readText: file => inspection.readFile(file, true),
      cache: inputMetadata,
    }).then(validation => {
      if (inspection.signal.aborted) return validation;
      inputsReady = validation.ready;
      renderInputValidation(validation);
      renderRecordRanges(validation.filesWithRecords);
      onChange();
      return validation;
    }).catch(error => {
      const validation = {
        ready: false,
        rows: [{
          name: "Selected inputs",
          detail: "Not ready",
          reason: error.message || String(error),
          ready: false,
        }],
      };
      if (!inspection.signal.aborted) {
        inputsReady = false;
        onChange();
        renderInputValidation(validation);
        emptyRecordRanges("Could not inspect selected inputs.");
      }
      return validation;
    }).finally(() => {
      if (inputInspection === inspection) {
        inputInspection = null;
        recordRangesElement.removeAttribute("aria-busy");
      }
      inspection.dispose();
    });
    onInputsChanged();
  }

  prefilterInput.addEventListener("change", syncPrefilterControls);
  noAlignInput.addEventListener("change", syncPrefilterControls);
  separateRecordsInput.addEventListener("change", onInputsChanged);
  filesInput.addEventListener("change", () => {
    setSelectedInputFiles([...filesInput.files]);
    handleFileSelection();
  });
  clearInputFilesButton.addEventListener("click", clear);
  let fileDragDepth = 0;
  fileDropZone.addEventListener("dragenter", event => {
    if (filesInput.disabled || ![...event.dataTransfer.types].includes("Files")) return;
    event.preventDefault();
    fileDragDepth += 1;
    fileDropZone.classList.add("is-dragging");
  });
  fileDropZone.addEventListener("dragover", event => {
    if (filesInput.disabled || ![...event.dataTransfer.types].includes("Files")) return;
    event.preventDefault();
    event.dataTransfer.dropEffect = "copy";
  });
  fileDropZone.addEventListener("dragleave", () => {
    fileDragDepth = Math.max(0, fileDragDepth - 1);
    if (fileDragDepth === 0) fileDropZone.classList.remove("is-dragging");
  });
  fileDropZone.addEventListener("drop", event => {
    event.preventDefault();
    fileDragDepth = 0;
    fileDropZone.classList.remove("is-dragging");
    if (filesInput.disabled) return;
    const files = [...event.dataTransfer.files];
    const unsupported = files.filter(file => !isSupportedInputFile(file));
    if (files.length === 0) return;
    if (unsupported.length > 0) {
      onError(
        `Unsupported file${unsupported.length === 1 ? "" : "s"}: ${unsupported.map(file => file.name).join(", ")}. Choose GenBank, GFF3/GFF, or FASTA files.`,
      );
      return;
    }
    setSelectedInputFiles(files);
    handleFileSelection();
  });
  geneFunctionsInput.addEventListener("change", syncAnnotationFiles);
  colourMapInput.addEventListener("change", syncAnnotationFiles);
  syncPrefilterControls();
  syncAnnotationFiles();

  function clear() {
    setSelectedInputFiles([], { replace: true });
    handleFileSelection();
  }

  function setDisabled(value) {
    disabled = value;
    for (const control of [filesInput, identityInput, noAlignInput, separateRecordsInput,
      keepOriginInput, recordRangesElement, geneFunctionsInput, clearInputFilesButton,
      container.querySelector("#run-tree-fields")]) control.disabled = disabled;
    fileInputStatusList.querySelectorAll(".remove-input-file").forEach(button => { button.disabled = disabled; });
    syncPrefilterControls();
    syncAnnotationFiles();
  }

  async function inspect() {
    let pending, validation;
    do {
      pending = inputInspectionPromise;
      validation = await pending;
    } while (pending !== inputInspectionPromise);
    return validation;
  }

  async function readOptions(workflow) {
    const identityCutoff = Number(identityInput.value) / 100;
    if (!identityInput.value.trim() || !Number.isFinite(identityCutoff) || identityCutoff < 0 || identityCutoff > 1) {
      throw new Error("Minimum protein identity must be a percentage between 0 and 100.");
    }
    const enabled = prefilterInput.checked && !noAlignInput.checked;
    const kmerSize = Number(kmerSizeInput.value);
    const minSharedKmers = Number(minSharedKmersInput.value);
    if (enabled && (!Number.isInteger(kmerSize) || kmerSize < 1 || kmerSize > 8)) {
      throw new Error("k-mer size must be a whole number between 1 and 8.");
    }
    if (enabled && (!Number.isInteger(minSharedKmers) || minSharedKmers < 1)) {
      throw new Error("Minimum shared k-mers must be a positive whole number.");
    }
    const settings = {
      enabled, noAlign: noAlignInput.checked, kmerSize, minSharedKmers, identityCutoff,
      genbank: { ranges: collectRecordRanges(), asSeparateClusters: separateRecordsInput.checked, setOrigin: !keepOriginInput.checked },
    };
    const useFileOrder = settings.noAlign || fileOrderInput.checked;
    const groupDefinitions = await readGroupDefinitions(workflow);
    return { settings, useFileOrder, groupDefinitions };
  }

  return {
    get ready() { return inputsReady; },
    get files() { return [...filesInput.files]; },
    setDisabled,
    inspect,
    readOptions,
    async treeInputs() {
      const validation = await inspect();
      return { ready: validation.ready, clusters: inputTreeClusters(validation.filesWithRecords || [], separateRecordsInput.checked) };
    },
    clear,
    focus() { filesInput.focus(); },
  };
}
