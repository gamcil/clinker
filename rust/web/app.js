import { displayedClusterOrder, similarityCsv, parseProject, downloadText } from "./project.mjs?v=20261015-compact-kmers-v1";
import { classifyInputFiles, fileStem, isGenBankFile, isGffFile, isFastaFile, isSupportedInputFile } from "./input-files.mjs?v=20261015-start-v35";
import { Workflow } from "./workflow.mjs?v=20261015-start-v35";
import {
  defineClinkerEditor,
  mountPlotSelectionToolbar,
} from "./clustermap-editor.mjs?v=20261014-clustermap-905a0f8-v16";

defineClinkerEditor();

const filesInput = document.querySelector("#files");
const fileDropZone = document.querySelector("#file-drop-zone");
const fileSummary = document.querySelector("#file-summary");
const fileInputStatus = document.querySelector("#file-input-status");
const fileInputStatusSummary = document.querySelector("#file-input-status-summary");
const fileInputStatusList = document.querySelector("#file-input-status-list");
const clearInputFilesButton = document.querySelector("#clear-input-files");
const identityInput = document.querySelector("#identity");
const prefilterInput = document.querySelector("#prefilter");
const noAlignInput = document.querySelector("#no-align");
const fileOrderInput = document.querySelector("#file-order");
const kmerSizeInput = document.querySelector("#kmer-size");
const minSharedKmersInput = document.querySelector("#min-shared-kmers");
const separateRecordsInput = document.querySelector("#separate-records");
const keepOriginInput = document.querySelector("#keep-origin");
const recordRangesElement = document.querySelector("#record-ranges");
const geneFunctionsInput = document.querySelector("#gene-functions");
const colourMapInput = document.querySelector("#colour-map");
const geneFunctionsSummary = document.querySelector("#gene-functions-summary");
const colourMapSummary = document.querySelector("#colour-map-summary");
const analyseButton = document.querySelector("#analyse");
const exampleButton = document.querySelector("#load-example");
const matrixButton = document.querySelector("#download-matrix");
const svgButton = document.querySelector("#download-svg");
const projectButton = document.querySelector("#download-project");
const loadClinkerFileButton = document.querySelector("#load-clinker-file");
const openProjectStartButton = document.querySelector("#open-project-start");
const plotActions = document.querySelector("#plot-actions");
const clinkerFileInput = document.querySelector("#clinker-file");
const plotEditor = document.querySelector("#plot-editor");
const editPlotButton = document.querySelector("#edit-plot");
const runPanel = document.querySelector("#run-panel");
const newRunButton = document.querySelector("#new-run");
const cancelRunButton = document.querySelector("#cancel-run");
const runFeedback = document.querySelector("#run-feedback");
const runStepsElement = document.querySelector("#run-steps");
const runSteps = new Map([...document.querySelectorAll("[data-run-step]")].map(step => [
  step.dataset.runStep,
  {
    element: step,
    result: step.querySelector(".step-result"),
    progress: step.querySelector("progress"),
  },
]));
const runStatus = document.querySelector("#run-status");
const plotStatus = document.querySelector("#plot-status");
const plot = d3.select("#plot");
const plotStage = document.querySelector(".plot-stage");
const toolbarMenus = [...document.querySelectorAll(".menu")];
const GENES_PER_TILE_SIDE = 40;
const PAIRS_PER_TILE = 1600;
const WORKER_VERSION = "20261015-compact-kmers-v1";
const WORKER_URL = new URL("./worker.js", import.meta.url);
WORKER_URL.searchParams.set("v", WORKER_VERSION);
const BASE_CHART_CONFIG = {
  link: { bestOnly: true },
  plot: { renderer: "webgpu" },
  legend: { columns: 4, position: "bottom" },
};
function createChart(config = {}) {
  return ClusterMap.ClusterMap().config(BASE_CHART_CONFIG).config(config);
}

// A chart keeps its camera and editable locus state. Each completed analysis
// is a new dataset, so give it a fresh instance rather than carrying state
// between positional IDs such as `cluster-0` and `locus-0-0`.
let chart = createChart();
plotEditor.chart = chart;
let disposeSelectionToolbar = null;
let latestSimilarity = null;
let latestAnalysisSession = null;
let analysisInProgress = false;
let activeWorkflow = null;
let recordInspectionVersion = 0;
let recordInspectionPromise = Promise.resolve();
let inputValidationPromise = Promise.resolve({ ready: false, rows: [] });
let inputsReady = false;
const selectedInputFiles = new Map();
const EXAMPLE_FILES = [
  "P. vexata CBS 129021.gbk",
  "A. versicolor CBS 583.65.gbk",
  "A. mulundensis DSM 5745.gbk",
  "A. burnettii MST-FP2249.gbk",
  "A. alliaceus CBS 536.65.gbk",
];

function setAnalysisInputsDisabled(disabled) {
  filesInput.disabled = disabled;
  identityInput.disabled = disabled;
  prefilterInput.disabled = disabled;
  noAlignInput.disabled = disabled;
  fileOrderInput.disabled = disabled;
  separateRecordsInput.disabled = disabled;
  keepOriginInput.disabled = disabled;
  recordRangesElement.disabled = disabled;
  geneFunctionsInput.disabled = disabled;
  colourMapInput.disabled = disabled || geneFunctionsInput.files.length === 0;
  analyseButton.disabled = disabled || !inputsReady;
  clearInputFilesButton.disabled = disabled;
  fileInputStatusList.querySelectorAll(".remove-input-file").forEach(button => {
    button.disabled = disabled;
  });
  exampleButton.disabled = disabled;
  loadClinkerFileButton.disabled = disabled;
  openProjectStartButton.disabled = disabled;
  newRunButton.disabled = disabled;
  cancelRunButton.disabled = false;
  cancelRunButton.hidden = !disabled && !chart.data();
  cancelRunButton.textContent = disabled ? "Cancel" : "Back to plot";
  syncPrefilterControls();
}

function syncPrefilterControls() {
  prefilterInput.disabled = noAlignInput.disabled || noAlignInput.checked;
  fileOrderInput.disabled = noAlignInput.disabled || noAlignInput.checked;
  const disabled = prefilterInput.disabled || !prefilterInput.checked;
  kmerSizeInput.disabled = disabled;
  minSharedKmersInput.disabled = disabled;
}

async function runWorkflow(task, failureMessage) {
  if (activeWorkflow) return;
  const workflow = new Workflow();
  const trigger = document.activeElement;
  openRunPanel();
  activeWorkflow = workflow;
  analysisInProgress = true;
  document.querySelector(".run-card").classList.add("is-busy");
  setAnalysisInputsDisabled(true);
  matrixButton.disabled = true;
  svgButton.disabled = true;
  projectButton.disabled = true;
  editPlotButton.disabled = true;
  cancelRunButton.focus();
  let completed = false;
  try {
    completed = await task(workflow) === true;
    workflow.throwIfCancelled();
  } catch (error) {
    const cancelled = workflow.signal.aborted;
    for (const step of runSteps.values()) {
      if (step.element.dataset.state !== "active") continue;
      step.element.dataset.state = cancelled ? "cancelled" : "error";
      step.progress.hidden = true;
    }
    setRunStatus(cancelled ? "Cancelled. You can change your settings and try again." : `${failureMessage}: ${error.message || String(error)}`, {
      tone: cancelled ? "" : "error",
      preserveSteps: true,
    });
    completed = false;
  } finally {
    workflow.dispose();
    activeWorkflow = null;
    analysisInProgress = false;
    document.querySelector(".run-card").classList.remove("is-busy");
    setAnalysisInputsDisabled(false);
    const hasPlot = Boolean(chart.data());
    matrixButton.disabled = !hasPlot || !latestSimilarity;
    svgButton.disabled = !hasPlot;
    projectButton.disabled = !hasPlot || !latestAnalysisSession;
    editPlotButton.disabled = !hasPlot;
    if (completed) {
      closeRunPanel();
      editPlotButton.focus();
    } else if (trigger instanceof HTMLElement && !trigger.disabled && trigger.getClientRects().length) {
      trigger.focus();
    }
  }
}

function parseTwoColumnCsv(text, description) {
  const rows = [];
  let row = [];
  let field = "";
  let quoted = false;
  for (let index = 0; index <= text.length; index += 1) {
    const character = index === text.length ? "\n" : text[index];
    if (quoted) {
      if (character === '"' && text[index + 1] === '"') {
        field += '"';
        index += 1;
      } else if (character === '"') {
        quoted = false;
      } else {
        field += character;
      }
    } else if (character === '"' && field.length === 0) {
      quoted = true;
    } else if (character === ",") {
      row.push(field.trim());
      field = "";
    } else if (character === "\n") {
      row.push(field.trim());
      field = "";
      if (row.some(value => value.length > 0)) {
        if (row.length !== 2 || row.some(value => !value)) {
          throw new Error(`${description} row ${rows.length + 1} must contain two non-empty columns.`);
        }
        rows.push(row);
      }
      row = [];
    } else if (character !== "\r") {
      field += character;
    }
  }
  if (quoted) throw new Error(`${description} contains an unterminated quoted field.`);
  return rows;
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
    remove.className = "remove-input-file quiet-action";
    remove.textContent = "Remove";
    remove.setAttribute("aria-label", `Remove ${row.name}`);
    remove.disabled = analysisInProgress;
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

function genBankRecords(text) {
  const records = [];
  const locusLine = /^LOCUS\s+(\S+)\s+(\d+)\s+bp(?:\s|$)/gim;
  for (const match of text.matchAll(locusLine)) {
    const length = Number(match[2]);
    if (Number.isSafeInteger(length) && length > 0) {
      records.push({ name: match[1], start: 1, end: length, length });
    }
  }
  return records;
}

function fastaRecords(text) {
  const records = [];
  let current = null;
  for (const line of text.split(/\r?\n/)) {
    if (line.startsWith(">")) {
      if (current) records.push(current);
      current = { name: line.slice(1).trim().split(/\s+/, 1)[0], start: 1, end: 0, length: 0 };
    } else if (current) {
      current.length += line.replace(/\s/g, "").length;
    }
  }
  if (current) records.push(current);
  return records
    .filter(record => record.name && record.length > 0)
    .map(record => ({ ...record, end: record.length }));
}

function gffRecords(text, companionText = "") {
  const records = [...text.matchAll(/^##sequence-region\s+(\S+)\s+(\d+)\s+(\d+)\s*$/gim)]
    .map(match => ({
      name: match[1],
      start: Number(match[2]),
      end: Number(match[3]),
      length: Number(match[3]) - Number(match[2]) + 1,
    }))
    .filter(record => Number.isSafeInteger(record.start)
      && Number.isSafeInteger(record.end)
      && Number.isSafeInteger(record.length)
      && record.start > 0
      && record.end >= record.start);
  if (records.length > 0) return records;
  const embedded = text.match(/^##FASTA\s*$([\s\S]*)/im)?.[1] || companionText;
  return fastaRecords(embedded);
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
  headings.className = "sequence-regions sequence-region-headings";
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
    table.className = "sequence-regions";
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

async function inspectSelectedFiles(files, version) {
  recordRangesElement.setAttribute("aria-busy", "true");
  emptyRecordRanges(files.length ? "Reading annotation records…" : "Choose input files to trim their sequences.");
  const fastaTexts = new Map(await Promise.all(files.filter(isFastaFile).map(async file => [
    fileStem(file.name),
    await file.text(),
  ])));
  const inspected = await Promise.all(files
    .filter(file => isGenBankFile(file) || isGffFile(file))
    .map(async file => {
      const text = await file.text();
      return {
        file,
        records: isGenBankFile(file)
          ? genBankRecords(text)
          : gffRecords(text, fastaTexts.get(fileStem(file.name)) || ""),
      };
    }));
  if (version !== recordInspectionVersion) return;
  renderRecordRanges(inspected);
  recordRangesElement.removeAttribute("aria-busy");
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
  const version = ++recordInspectionVersion;
  const files = [...filesInput.files];
  inputsReady = false;
  analyseButton.disabled = true;
  fileInputStatusList.replaceChildren();
  fileInputStatus.hidden = files.length === 0;
  if (files.length > 0) fileInputStatusSummary.textContent = "Checking selected inputs…";
  inputValidationPromise = classifyInputFiles(files)
    .then(validation => {
      if (version !== recordInspectionVersion) return validation;
      inputsReady = validation.ready;
      renderInputValidation(validation);
      analyseButton.disabled = analysisInProgress || !inputsReady;
      return validation;
    })
    .catch(error => {
      const validation = {
        ready: false,
        rows: [{
          name: "Selected inputs",
          detail: "Not ready",
          reason: error.message || String(error),
          ready: false,
        }],
      };
      if (version === recordInspectionVersion) {
        inputsReady = false;
        renderInputValidation(validation);
        analyseButton.disabled = true;
      }
      return validation;
    });
  recordInspectionPromise = inspectSelectedFiles(files, version)
    .catch(error => {
      if (version !== recordInspectionVersion) return;
      emptyRecordRanges(`Could not inspect record names: ${error.message || String(error)}`);
      recordRangesElement.removeAttribute("aria-busy");
    });
}

prefilterInput.addEventListener("change", syncPrefilterControls);
noAlignInput.addEventListener("change", syncPrefilterControls);
filesInput.addEventListener("change", () => {
  setSelectedInputFiles([...filesInput.files]);
  handleFileSelection();
});
clearInputFilesButton.addEventListener("click", () => {
  setSelectedInputFiles([], { replace: true });
  handleFileSelection();
});
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
    setRunStatus(
      `Unsupported file${unsupported.length === 1 ? "" : "s"}: ${unsupported.map(file => file.name).join(", ")}. Choose GenBank, GFF3/GFF, or FASTA files.`,
      { tone: "error" },
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

function resetRunSteps() {
  for (const step of runSteps.values()) {
    step.element.dataset.state = "waiting";
    step.result.textContent = "Waiting";
    step.progress.hidden = true;
    step.progress.value = 0;
  }
}

function showRunSteps({ includePrefilter = false } = {}) {
  resetRunSteps();
  runSteps.get("prefilter").element.hidden = !includePrefilter;
  let visibleIndex = 0;
  for (const step of runSteps.values()) {
    if (step.element.hidden) continue;
    visibleIndex += 1;
    step.element.querySelector(".step-marker").textContent = visibleIndex;
  }
  runFeedback.hidden = false;
  runStepsElement.hidden = false;
  runStatus.hidden = true;
  runStatus.textContent = "";
  delete runStatus.dataset.tone;
}

function setRunStep(name, state, result, { progress = null } = {}) {
  const step = runSteps.get(name);
  step.element.dataset.state = state;
  step.result.textContent = result;
  step.progress.hidden = progress === null;
  if (progress === "indeterminate") step.progress.removeAttribute("value");
  else if (progress !== null) step.progress.value = progress;
}

function setRunStatus(message, { tone = "", preserveSteps = false } = {}) {
  runFeedback.hidden = !message;
  if (!preserveSteps) runStepsElement.hidden = true;
  runStatus.hidden = !message;
  runStatus.textContent = message;
  if (tone) runStatus.dataset.tone = tone;
  else delete runStatus.dataset.tone;
}

function setPlotStatus(message, { tone = "" } = {}) {
  plotStatus.hidden = !message;
  plotStatus.textContent = message;
  if (tone) plotStatus.dataset.tone = tone;
  else delete plotStatus.dataset.tone;
}

function setContextStatus(message, options = {}) {
  if (runPanel.hidden) setPlotStatus(message, options);
  else setRunStatus(message, options);
}

for (const menu of toolbarMenus) {
  menu.addEventListener("toggle", () => {
    if (!menu.open) return;
    toolbarMenus.forEach(other => {
      if (other !== menu) other.open = false;
    });
  });
  menu.addEventListener("click", event => {
    if (event.target.closest("button:not(:disabled)")) menu.open = false;
  });
}

document.addEventListener("pointerdown", event => {
  toolbarMenus.forEach(menu => {
    if (!menu.contains(event.target)) menu.open = false;
  });
});

document.addEventListener("keydown", event => {
  if (event.key !== "Escape") return;
  toolbarMenus.forEach(menu => { menu.open = false; });
  if (!analysisInProgress && !runPanel.hidden && chart.data()) closeRunPanel();
});

function closeRunPanel() {
  if (analysisInProgress || !chart.data()) return;
  runPanel.hidden = true;
  plotActions.hidden = false;
  newRunButton.hidden = false;
  newRunButton.setAttribute("aria-expanded", "false");
  cancelRunButton.hidden = true;
}

function openRunPanel() {
  runPanel.hidden = false;
  newRunButton.setAttribute("aria-expanded", "true");
  cancelRunButton.hidden = !chart.data();
  cancelRunButton.textContent = "Back to plot";
  setRunStatus("");
  closePlotEditor();
}

newRunButton.addEventListener("click", () => {
  if (runPanel.hidden) openRunPanel();
  else closeRunPanel();
});
cancelRunButton.addEventListener("click", () => {
  if (activeWorkflow) activeWorkflow.cancel();
  else closeRunPanel();
});

function closePlotEditor() {
  plotEditor.hidden = true;
  editPlotButton.setAttribute("aria-expanded", "false");
}

function showPlotEditor() {
  if (!runPanel.hidden) closeRunPanel();
  plotEditor.hidden = false;
  editPlotButton.setAttribute("aria-expanded", "true");
}

editPlotButton.addEventListener("click", () => {
  const isOpen = !plotEditor.hidden;
  if (isOpen) closePlotEditor();
  else showPlotEditor();
});

function workerCount(taskCount) {
  // Each worker owns a separate Wasm instance and its alignment buffers. A
  // modest cap avoids making many large DP matrices resident at once.
  return Math.min(taskCount, Math.max(1, Math.min(navigator.hardwareConcurrency || 2, 4)));
}

function parseInputFiles(
  files,
  prefilter,
  { onProgress, onPrefilterStart, onPrefilterProgress },
  workflow,
) {
  const worker = workflow.worker(WORKER_URL);
  // These buffers belong to this run; transfer them once without another copy.
  const taskFiles = files;

  return new Promise((resolve, reject) => {
    let nextFileIndex = 0;
    let settled = false;
    const startupTimeout = setTimeout(() => {
      finish(reject, new Error("The analysis worker did not start. Reload the page and try again."));
    }, 30000);
    const finish = (callback, value) => {
      if (settled) return;
      settled = true;
      clearTimeout(startupTimeout);
      workflow.signal.removeEventListener("abort", cancel);
      workflow.releaseWorker(worker);
      callback(value);
    };
    const cancel = () => finish(reject, workflow.signal.reason);
    workflow.signal.addEventListener("abort", cancel, { once: true });
    const sendNextFile = () => {
      if (settled) return;
      if (nextFileIndex === taskFiles.length) {
        if (prefilter.enabled) onPrefilterStart();
        worker.postMessage({ type: "parse-finish" });
        return;
      }
      const fileIndex = nextFileIndex;
      const file = taskFiles[nextFileIndex++];
      worker.postMessage(
        { type: "parse-file", file, fileIndex },
        [file.bytes.buffer],
      );
    };
    const continueAfterPaint = () => {
      if (document.hidden) {
        setTimeout(sendNextFile, 0);
      } else {
        requestAnimationFrame(() => setTimeout(sendNextFile, 0));
      }
    };
    worker.onmessage = ({ data }) => {
      if (workflow.signal.aborted) return;
      if (data.type === "error") finish(reject, new Error(data.message));
      else if (data.type === "parse-ready") {
        clearTimeout(startupTimeout);
        sendNextFile();
      }
      else if (data.type === "file-parsed") {
        onProgress(data.fileIndex + 1, taskFiles.length);
        continueAfterPaint();
      } else if (data.type === "filter-progress") {
        onPrefilterProgress(data.completed, data.total);
      } else if (data.type === "parsed") {
        finish(resolve, { ...data.parsed, wasmModule: data.wasmModule });
      }
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage({ type: "parse-start", prefilter });
  });
}

function* candidateTiles(proteins, candidatePairs) {
  for (let start = 0; start < candidatePairs.length; start += PAIRS_PER_TILE) {
    const proteinIndexes = new Map();
    const tileProteins = [];
    const localIndex = index => {
      if (!proteinIndexes.has(index)) {
        proteinIndexes.set(index, tileProteins.length);
        tileProteins.push(proteins[index]);
      }
      return proteinIndexes.get(index);
    };
    const pairs = candidatePairs.slice(start, start + PAIRS_PER_TILE).map(pair => ({
      queryIndex: localIndex(pair.queryIndex),
      targetIndex: localIndex(pair.targetIndex),
    }));
    yield { type: "analyse-pairs", proteins: tileProteins, pairs, pairCount: pairs.length };
  }
}

function proteinsByCluster(proteins, clusterCount) {
  const clusters = Array.from({ length: clusterCount }, () => []);
  proteins.forEach(protein => clusters[protein.cluster].push(protein));
  return clusters;
}

function tileCount(clusters) {
  let count = 0;
  for (let queryCluster = 0; queryCluster < clusters.length; queryCluster += 1) {
    for (let targetCluster = queryCluster + 1; targetCluster < clusters.length; targetCluster += 1) {
      count += Math.ceil(clusters[queryCluster].length / GENES_PER_TILE_SIDE)
        * Math.ceil(clusters[targetCluster].length / GENES_PER_TILE_SIDE);
    }
  }
  return count;
}

function pairCount(clusters) {
  let count = 0;
  for (let queryCluster = 0; queryCluster < clusters.length; queryCluster += 1) {
    for (let targetCluster = queryCluster + 1; targetCluster < clusters.length; targetCluster += 1) {
      count += clusters[queryCluster].length * clusters[targetCluster].length;
    }
  }
  return count;
}

function progressPercent(completed, total) {
  if (total === 0) return 100;
  return Math.min(100, Math.floor((completed / total) * 100));
}

function prefilterResult(aligned, possible) {
  const skipped = Math.max(0, possible - aligned);
  const percentage = possible === 0 ? 0 : (skipped / possible) * 100;
  const formattedPercentage = percentage < 10 ? percentage.toFixed(1) : Math.round(percentage);
  return `Filtered out ${skipped.toLocaleString()} of ${possible.toLocaleString()} possible alignments (${formattedPercentage}%)`;
}

function* alignmentTiles(clusters) {
  for (let queryCluster = 0; queryCluster < clusters.length; queryCluster += 1) {
    for (let targetCluster = queryCluster + 1; targetCluster < clusters.length; targetCluster += 1) {
      const queryProteins = clusters[queryCluster];
      const targetProteins = clusters[targetCluster];
      for (let queryStart = 0; queryStart < queryProteins.length; queryStart += GENES_PER_TILE_SIDE) {
        for (let targetStart = 0; targetStart < targetProteins.length; targetStart += GENES_PER_TILE_SIDE) {
          yield {
            type: "analyse-tile",
            query: queryProteins.slice(queryStart, queryStart + GENES_PER_TILE_SIDE),
            target: targetProteins.slice(targetStart, targetStart + GENES_PER_TILE_SIDE),
            pairCount: Math.min(GENES_PER_TILE_SIDE, queryProteins.length - queryStart)
              * Math.min(GENES_PER_TILE_SIDE, targetProteins.length - targetStart),
          };
        }
      }
    }
  }
}

function plotLink(link) {
  return {
    // Keep positional references too: the grouping worker needs them to build
    // Rust-side union-find components, while clustermap.js uses only `uid`.
    query: {
      ...link.query,
      uid: `gene-${link.query.cluster}-${link.query.locus}-${link.query.gene}`,
    },
    target: {
      ...link.target,
      uid: `gene-${link.target.cluster}-${link.target.locus}-${link.target.gene}`,
    },
    identity: link.identity,
    similarity: link.similarity,
  };
}

function bestLinkCount(plotData) {
  // This mirrors the library's best-only rule for the status summary. The
  // renderer still receives every link so the appearance control can switch
  // between the reduced and complete link sets at runtime.
  const clusterForGene = new Map();
  plotData.clusters.forEach(cluster => {
    cluster.loci.forEach(locus => {
      locus.genes.forEach(gene => clusterForGene.set(gene.uid, cluster.uid));
    });
  });
  const pairKey = (left, right) => [left, right].sort().join("\u0000");
  const selectedByPair = new Map();

  for (const link of [...plotData.links].sort((left, right) => right.identity - left.identity)) {
    const queryCluster = clusterForGene.get(link.query.uid);
    const targetCluster = clusterForGene.get(link.target.uid);
    if (!queryCluster || !targetCluster) continue;

    const key = pairKey(queryCluster, targetCluster);
    const selected = selectedByPair.get(key) || [];
    const superseded = selected.some(candidate =>
      link.identity < candidate.identity
      && (candidate.query.uid === link.query.uid
        || candidate.query.uid === link.target.uid
        || candidate.target.uid === link.query.uid
        || candidate.target.uid === link.target.uid),
    );
    if (!superseded) selected.push(link);
    selectedByPair.set(key, selected);
  }

  return [...selectedByPair.values()].reduce((count, links) => count + links.length, 0);
}

matrixButton.addEventListener("click", () => {
  if (!latestSimilarity) return;
  try {
    downloadText("clinker-cluster-similarity.csv", similarityCsv(latestSimilarity, displayedClusterOrder({ data: chart.data(), state: chart.state() })), "text/csv;charset=utf-8");
  } catch (error) {
    setContextStatus(`CSV export failed: ${error.message || String(error)}`, { tone: "error" });
  }
});

svgButton.addEventListener("click", () => {
  try {
    downloadText("clinker.svg", chart.exportSvg(), "image/svg+xml;charset=utf-8");
  } catch (error) {
    setContextStatus(`SVG export failed: ${error.message || String(error)}`, { tone: "error" });
  }
});

function replaceChart(config = chart.config()) {
  const previousChart = chart;
  disposeSelectionToolbar?.();
  disposeSelectionToolbar = null;
  chart = createChart(config);
  previousChart.destroy();
}

function mountSelectionToolbar() {
  disposeSelectionToolbar?.();
  disposeSelectionToolbar = mountPlotSelectionToolbar(chart, plotStage);
}

function mountImportedPlot(data, { config, state = null }) {
  replaceChart(config);
  plot.datum(data).call(chart);
  if (state) chart.state(state);
  mountSelectionToolbar();
  plotEditor.chart = chart;
  latestSimilarity = null;
  latestAnalysisSession = null;
  matrixButton.disabled = true;
  svgButton.disabled = false;
  projectButton.disabled = true;
  editPlotButton.disabled = false;
  closePlotEditor();
  closeRunPanel();
}

projectButton.addEventListener("click", () => {
  if (!latestAnalysisSession) return;
  const project = { ...chart.project(), version: 2, analysis: latestAnalysisSession };
  downloadText("clinker-project.json", `${JSON.stringify(project, null, 2)}\n`, "application/json;charset=utf-8");
});

function openAnalysisSession(session, workflow) {
  const worker = workflow.worker(WORKER_URL);
  return new Promise((resolve, reject) => {
    const finish = (callback, value) => {
      workflow.signal.removeEventListener("abort", cancel);
      workflow.releaseWorker(worker);
      callback(value);
    };
    const cancel = () => finish(reject, workflow.signal.reason);
    workflow.signal.addEventListener("abort", cancel, { once: true });
    worker.onmessage = ({ data }) => {
      if (workflow.signal.aborted) return;
      if (data.type === "error") finish(reject, new Error(data.message));
      else if (data.type === "post-process") finish(resolve, data.result);
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage({ type: "open-session", session, useFileOrder: false });
  });
}

function retainAnalysisResult(result) {
  latestSimilarity = {
    clusterNames: result.clusterNames,
    similarityMatrix: result.similarityMatrix,
  };
  latestAnalysisSession = result.analysisSession;
  matrixButton.disabled = false;
  projectButton.disabled = false;
}

loadClinkerFileButton.addEventListener("click", () => clinkerFileInput.click());
openProjectStartButton.addEventListener("click", () => clinkerFileInput.click());
clinkerFileInput.addEventListener("change", () => {
  const file = clinkerFileInput.files[0];
  clinkerFileInput.value = "";
  if (!file) return;
  runWorkflow(async workflow => {
    setRunStatus("Opening project…");
    const saved = parseProject(await workflow.readFile(file, true));
    workflow.throwIfCancelled();
    const analysisResult = await openAnalysisSession(saved.analysis, workflow);
    workflow.throwIfCancelled();
    mountImportedPlot(saved.data, { config: saved.config, state: saved.state });
    retainAnalysisResult(analysisResult);
    setPlotStatus("Loaded clinker project with its analysis data.");
    return true;
  }, "Could not open project");
});

function postProcess(
  layout,
  proteins,
  links,
  groupDefinitions,
  identity,
  useFileOrder,
  wasmModule,
  onProgress,
  workflow,
) {
  const worker = workflow.worker(WORKER_URL);
  return new Promise((resolve, reject) => {
    const finish = (callback, value) => {
      workflow.signal.removeEventListener("abort", cancel);
      workflow.releaseWorker(worker);
      callback(value);
    };
    const cancel = () => finish(reject, workflow.signal.reason);
    workflow.signal.addEventListener("abort", cancel, { once: true });
    worker.onmessage = ({ data }) => {
      if (workflow.signal.aborted) return;
      if (data.type === "error") finish(reject, new Error(data.message));
      else if (data.type === "post-process-progress") onProgress(data);
      else if (data.type === "post-process") finish(resolve, data.result);
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage({
      type: "post-process",
      layout,
      proteins,
      useFileOrder,
      wasmModule,
      identity,
      reportProgress: true,
      links: links.map(link => ({
        query: link.query,
        target: link.target,
        identity: link.identity,
        similarity: link.similarity,
      })),
      groupDefinitions,
    });
  });
}

function analyseTilesInWorkerPool(
  proteins,
  clusterCount,
  identity,
  candidatePairs,
  wasmModule,
  onProgress,
  workflow,
) {
  const clusters = proteinsByCluster(proteins, clusterCount);
  const totalTiles = candidatePairs
    ? Math.ceil(candidatePairs.length / PAIRS_PER_TILE)
    : tileCount(clusters);
  if (totalTiles === 0) return Promise.resolve({ links: [], alignmentCount: 0 });
  const tiles = candidatePairs ? candidateTiles(proteins, candidatePairs) : alignmentTiles(clusters);
  const totalPairs = candidatePairs ? candidatePairs.length : pairCount(clusters);
  const workers = Array.from(
    { length: workerCount(totalTiles) },
    () => workflow.worker(WORKER_URL),
  );
  const linksByTask = new Array(totalTiles);
  const alignmentCountsByTask = new Array(totalTiles).fill(0);
  const processedPairsByTask = new Array(totalTiles).fill(0);
  let nextTask = 0;
  let completed = 0;
  let processedPairs = 0;
  let settled = false;

  return new Promise((resolve, reject) => {
    const finish = error => {
      if (settled) return;
      settled = true;
      workflow.signal.removeEventListener("abort", cancel);
      workers.forEach(worker => workflow.releaseWorker(worker));
      if (error) reject(error);
      else {
        const links = linksByTask.flat();
        links.forEach((link, index) => { link.uid = `link-${index}`; });
        resolve({
          links,
          alignmentCount: alignmentCountsByTask.reduce((sum, count) => sum + count, 0),
        });
      }
    };

    const cancel = () => finish(workflow.signal.reason);
    workflow.signal.addEventListener("abort", cancel, { once: true });
    const reportTaskProgress = (taskIndex, count) => {
      const next = Math.max(processedPairsByTask[taskIndex], count);
      processedPairs += next - processedPairsByTask[taskIndex];
      processedPairsByTask[taskIndex] = next;
      onProgress(processedPairs, totalPairs);
    };
    const dispatch = worker => {
      if (nextTask === totalTiles) return;
      const taskIndex = nextTask++;
      const tile = tiles.next().value;
      processedPairsByTask[taskIndex] = 0;
      worker.postMessage({ ...tile, identity, taskIndex, wasmModule });
    };

    workers.forEach(worker => {
      worker.onmessage = ({ data }) => {
        if (settled) return;
        if (data.type === "error") {
          finish(new Error(data.message));
          return;
        }

        if (data.type === "tile-progress") {
          reportTaskProgress(data.taskIndex, data.processedPairs);
          return;
        }

        if (data.type !== "tile-result") return;
        linksByTask[data.taskIndex] = data.links.map(plotLink);
        alignmentCountsByTask[data.taskIndex] = data.alignmentCount || 0;
        reportTaskProgress(
          data.taskIndex,
          data.pairCount || processedPairsByTask[data.taskIndex],
        );
        completed += 1;
        if (completed === totalTiles) finish();
        else dispatch(worker);
      };
      worker.onerror = event => finish(new Error(event.message || "Worker failed"));
      dispatch(worker);
    });
  });
}

async function analyseFiles(files, workflow) {
  const identity = Number(identityInput.value) / 100;
  if (!identityInput.value.trim() || !Number.isFinite(identity) || identity < 0 || identity > 1) {
    setRunStatus("Minimum protein identity must be a percentage between 0 and 100.", { tone: "error" });
    return;
  }
  const kmerSize = Number(kmerSizeInput.value);
  const minSharedKmers = Number(minSharedKmersInput.value);
  if (!noAlignInput.checked && prefilterInput.checked
      && (!Number.isInteger(kmerSize) || kmerSize < 1 || kmerSize > 8)) {
    setRunStatus("k-mer size must be a whole number between 1 and 8.", { tone: "error" });
    return;
  }
  if (!noAlignInput.checked && prefilterInput.checked
      && (!Number.isInteger(minSharedKmers) || minSharedKmers < 1)) {
    setRunStatus("Minimum shared k-mers must be a positive whole number.", { tone: "error" });
    return;
  }
  const ranges = collectRecordRanges();
  const groupDefinitions = await readGroupDefinitions(workflow);
  workflow.throwIfCancelled();
  const prefilter = {
    enabled: prefilterInput.checked && !noAlignInput.checked,
    kmerSize,
    minSharedKmers,
    identityCutoff: identity,
    genbank: {
      ranges,
      asSeparateClusters: separateRecordsInput.checked,
      setOrigin: !keepOriginInput.checked,
    },
  };
  showRunSteps({ includePrefilter: prefilter.enabled });
  setRunStep("parse", "active", `Parsing 0 of ${files.length} files…`, { progress: 0 });
  const parsed = await parseInputFiles(files, prefilter, {
    onProgress: (processed, total) => {
      setRunStep("parse", "active", `Parsed ${processed} of ${total} files…`, {
        progress: progressPercent(processed, total),
      });
    },
    onPrefilterStart: () => {
      setRunStep(
        "parse",
        "complete",
        `Parsed ${files.length} ${files.length === 1 ? "file" : "files"}`,
      );
      setRunStep("prefilter", "active", "Selecting candidate alignments…", {
        progress: 0,
      });
    },
    onPrefilterProgress: (processed, total) => {
      const progress = progressPercent(processed, total);
      setRunStep(
        "prefilter",
        "active",
        `Examined ${processed.toLocaleString()} of ${total.toLocaleString()} possible alignments…`,
        { progress },
      );
    },
  }, workflow);
  workflow.throwIfCancelled();
  setRunStep(
    "parse",
    "complete",
    `Parsed ${files.length} ${files.length === 1 ? "file" : "files"} · ${parsed.proteins.length.toLocaleString()} genes`,
  );
  const proteinsByInputCluster = proteinsByCluster(
    parsed.proteins,
    parsed.layout.clusters.length,
  );
  const possibleAlignments = pairCount(proteinsByInputCluster);
  const plannedAlignments = parsed.candidatePairs?.length ?? possibleAlignments;
  if (prefilter.enabled) {
    setRunStep(
      "prefilter",
      "complete",
      prefilterResult(plannedAlignments, possibleAlignments),
    );
  }
  setRunStep("align", "active", "Aligning genes… 0%", { progress: 0 });
  const alignmentResult = noAlignInput.checked
    ? { links: [], alignmentCount: 0 }
    : await analyseTilesInWorkerPool(
      parsed.proteins,
      parsed.layout.clusters.length,
      identity,
      parsed.candidatePairs,
      parsed.wasmModule,
      (completedPairs, total) => {
        const progress = progressPercent(completedPairs, total);
        setRunStep("align", "active", `Aligning genes… ${progress}%`, { progress });
      },
      workflow,
    );
  const { links, alignmentCount } = alignmentResult;
  setRunStep(
    "align",
    "complete",
    noAlignInput.checked
      ? "Skipped protein alignment"
      : `Ran ${alignmentCount.toLocaleString()} alignments · retained ${links.length.toLocaleString()} links`,
  );
  const useFileOrder = noAlignInput.checked || fileOrderInput.checked;
  setRunStep("order", "active", "Scoring cluster pairs…", { progress: 0 });
  const result = await postProcess(
    parsed.layout,
    parsed.proteins,
    links,
    groupDefinitions,
    identity,
    useFileOrder,
    parsed.wasmModule,
    ({ stage, completed, total }) => {
      if (stage === "similarity") {
        setRunStep(
          "order",
          "active",
          `Scored ${completed.toLocaleString()} of ${total.toLocaleString()} cluster pairs…`,
          { progress: progressPercent(completed, total) },
        );
      } else if (stage === "ordering") {
        if (completed < total) {
          setRunStep(
            "order",
            "active",
            useFileOrder ? "Applying input file order…" : "Ordering clusters by similarity…",
          );
        } else {
          setRunStep(
            "order",
            "complete",
            useFileOrder
              ? `Kept ${total.toLocaleString()} clusters in input order`
              : `Ordered ${total.toLocaleString()} clusters by similarity`,
          );
        }
      } else if (stage === "groups") {
        const linkTotal = total / 2;
        if (completed >= total) {
          setRunStep(
            "groups",
            "complete",
            `Grouped ${linkTotal.toLocaleString()} links`,
          );
          return;
        }
        const indexing = completed <= linkTotal;
        const phaseCompleted = indexing ? completed : completed - linkTotal;
        setRunStep(
          "groups",
          "active",
          indexing
            ? `Indexed ${phaseCompleted.toLocaleString()} of ${linkTotal.toLocaleString()} links…`
            : `Grouped ${phaseCompleted.toLocaleString()} of ${linkTotal.toLocaleString()} links…`,
          { progress: progressPercent(completed, total) },
        );
      } else if (stage === "layout") {
        setRunStep(
          "layout",
          "active",
          `Arranged ${completed.toLocaleString()} of ${total.toLocaleString()} clusters…`,
          { progress: progressPercent(completed, total) },
        );
      }
    },
    workflow,
  );
  workflow.throwIfCancelled();
  const plotData = result.plotData;
  latestSimilarity = {
    clusterNames: result.clusterNames,
    similarityMatrix: result.similarityMatrix,
  };
  latestAnalysisSession = result.analysisSession;
  matrixButton.disabled = false;
  const displayedBestLinks = bestLinkCount(plotData);
  const summary = `${plotData.clusters.length} clusters; ${plotData.links.length} retained links; ${displayedBestLinks} displayed best links.`;
  setRunStep(
    "groups",
    "complete",
    `Built ${plotData.groups.length.toLocaleString()} homology groups from ${links.length.toLocaleString()} links`,
  );
  const locusCount = plotData.clusters.reduce(
    (count, cluster) => count + cluster.loci.length,
    0,
  );
  setRunStep(
    "layout",
    "complete",
    `Arranged ${plotData.clusters.length.toLocaleString()} clusters · ${locusCount.toLocaleString()} loci`,
  );

  replaceChart();
  plot.datum(plotData).call(chart);
  mountSelectionToolbar();
  plotEditor.chart = chart;
  svgButton.disabled = false;
  projectButton.disabled = false;
  editPlotButton.disabled = false;
  setPlotStatus(summary);
  return true;
}

analyseButton.addEventListener("click", () => {
  if (filesInput.files.length === 0) return;
  const selected = [...filesInput.files];
  runWorkflow(async workflow => {
    setRunStatus("Preparing selected files…");
    const validation = await workflow.wait(inputValidationPromise);
    if (!validation.ready) throw new Error("Resolve the selected input-file problems before starting the analysis.");
    await workflow.wait(recordInspectionPromise);
    const files = [];
    for (const file of selected) {
      files.push({ name: file.name, bytes: new Uint8Array(await workflow.readFile(file)) });
    }
    return analyseFiles(files, workflow);
  }, "Analysis failed");
});

exampleButton.addEventListener("click", () => {
  runWorkflow(async workflow => {
    setRunStatus("Loading example files…");
    const files = await workflow.wait(Promise.all(EXAMPLE_FILES.map(async name => {
      const response = await fetch(`examples/${encodeURIComponent(name)}`, { signal: workflow.signal });
      if (!response.ok) throw new Error(`Could not load ${name} (${response.status})`);
      return { name, bytes: new Uint8Array(await response.arrayBuffer()) };
    })));
    workflow.throwIfCancelled();
    return analyseFiles(files, workflow);
  }, "Example analysis failed");
});
