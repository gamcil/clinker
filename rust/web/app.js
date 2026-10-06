import {
  defineClinkerEditor,
  mountPlotSelectionToolbar,
} from "./clustermap-editor.mjs?v=20261002-bundled-d3-v6";

defineClinkerEditor();

const filesInput = document.querySelector("#files");
const fileSummary = document.querySelector("#file-summary");
const identityInput = document.querySelector("#identity");
const prefilterInput = document.querySelector("#prefilter");
const clusterOrderInput = document.querySelector("#cluster-order");
const kmerSizeInput = document.querySelector("#kmer-size");
const minSharedKmersInput = document.querySelector("#min-shared-kmers");
const analyseButton = document.querySelector("#analyse");
const exampleButton = document.querySelector("#load-example");
const matrixButton = document.querySelector("#download-matrix");
const svgButton = document.querySelector("#download-svg");
const dataButton = document.querySelector("#download-data");
const loadDataButton = document.querySelector("#load-data");
const dataFileInput = document.querySelector("#data-file");
const projectButton = document.querySelector("#download-project");
const loadProjectButton = document.querySelector("#load-project");
const projectFileInput = document.querySelector("#project-file");
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
const WORKER_VERSION = "20261006-wasm-init-v11";
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
let analysisInProgress = false;
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
  clusterOrderInput.disabled = disabled;
  analyseButton.disabled = disabled;
  exampleButton.disabled = disabled;
  loadDataButton.disabled = disabled;
  loadProjectButton.disabled = disabled;
  newRunButton.disabled = disabled;
  cancelRunButton.disabled = disabled;
  syncPrefilterControls();
}

function syncPrefilterControls() {
  const disabled = prefilterInput.disabled || !prefilterInput.checked;
  kmerSizeInput.disabled = disabled;
  minSharedKmersInput.disabled = disabled;
}

function syncFileSummary() {
  const count = filesInput.files.length;
  fileSummary.textContent = count === 0
    ? "No files selected"
    : `${count} ${count === 1 ? "file" : "files"} selected`;
}

prefilterInput.addEventListener("change", syncPrefilterControls);
filesInput.addEventListener("change", syncFileSummary);
syncPrefilterControls();

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
  newRunButton.hidden = false;
  newRunButton.setAttribute("aria-expanded", "false");
  cancelRunButton.hidden = true;
}

function openRunPanel() {
  runPanel.hidden = false;
  newRunButton.setAttribute("aria-expanded", "true");
  cancelRunButton.hidden = !chart.data();
  setRunStatus("");
  closePlotEditor();
}

newRunButton.addEventListener("click", () => {
  if (runPanel.hidden) openRunPanel();
  else closeRunPanel();
});
cancelRunButton.addEventListener("click", closeRunPanel);

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

function copyFiles(files, fileIndexes) {
  // A transferred ArrayBuffer can only be sent once. Each cluster appears in
  // several comparisons, so make a task-local copy before transferring it.
  return fileIndexes.map(index => ({
    name: files[index].name,
    bytes: files[index].bytes.slice(),
  }));
}

function parseInputFiles(
  files,
  prefilter,
  { onProgress, onPrefilterStart, onPrefilterProgress },
) {
  const worker = new Worker(WORKER_URL, { type: "module" });
  const taskFiles = copyFiles(files, files.map((_, index) => index));

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
      worker.terminate();
      callback(value);
    };
    const sendNextFile = () => {
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
      } else if (data.type === "parsed") finish(resolve, data.parsed);
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

function csvCell(value) {
  const text = String(value);
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

function downloadSimilarityCsv({ clusterNames, clusterOrder, similarityMatrix }) {
  const plotPosition = new Map(clusterOrder.map((cluster, index) => [cluster, index + 1]));
  const rows = [[
    "query_plot_order",
    "query_cluster",
    "query_index",
    "target_plot_order",
    "target_cluster",
    "target_index",
    "query_to_target_coverage",
    "target_to_query_coverage",
    "containment_similarity",
    "distance",
  ]];
  for (let query = 0; query < clusterNames.length; query += 1) {
    for (let target = query + 1; target < clusterNames.length; target += 1) {
      const pair = similarityMatrix[query][target];
      rows.push([
        plotPosition.get(query),
        clusterNames[query],
        query,
        plotPosition.get(target),
        clusterNames[target],
        target,
        pair.queryCoverage,
        pair.targetCoverage,
        pair.similarity,
        1 - pair.similarity,
      ]);
    }
  }
  const csv = `${rows.map(row => row.map(csvCell).join(",")).join("\n")}\n`;
  const url = URL.createObjectURL(new Blob([csv], { type: "text/csv" }));
  const download = document.createElement("a");
  download.href = url;
  download.download = "clinker-cluster-similarity.csv";
  document.body.append(download);
  download.click();
  download.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

matrixButton.addEventListener("click", () => {
  if (latestSimilarity) downloadSimilarityCsv(latestSimilarity);
});

function downloadSvg() {
  const svg = chart.exportSvg();
  const url = URL.createObjectURL(new Blob([svg], { type: "image/svg+xml;charset=utf-8" }));
  const download = document.createElement("a");
  download.href = url;
  download.download = "clinker.svg";
  document.body.append(download);
  download.click();
  download.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}

svgButton.addEventListener("click", () => {
  try {
    downloadSvg();
  } catch (error) {
    setContextStatus(`SVG export failed: ${error.message || String(error)}`, { tone: "error" });
  }
});

dataButton.addEventListener("click", () => {
  const data = chart.data();
  if (!data) return;
  const url = URL.createObjectURL(new Blob([`${JSON.stringify(data, null, 2)}\n`], {
    type: "application/json;charset=utf-8",
  }));
  const download = document.createElement("a");
  download.href = url;
  download.download = "clinker-plot-data.json";
  document.body.append(download);
  download.click();
  download.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
});

function validatePlotData(data) {
  if (!data || typeof data !== "object" || Array.isArray(data)) {
    throw new Error("the file must contain a plot-data object");
  }
  if (!Array.isArray(data.clusters) || !Array.isArray(data.links) || !Array.isArray(data.groups)) {
    throw new Error("plot data must contain clusters, links, and groups arrays");
  }
  const genes = new Set();
  for (const cluster of data.clusters) {
    if (!cluster || !Array.isArray(cluster.loci)) throw new Error("every cluster needs a loci array");
    for (const locus of cluster.loci) {
      if (!locus || !Array.isArray(locus.genes)) throw new Error("every locus needs a genes array");
      for (const gene of locus.genes) {
        if (gene?.uid === undefined || genes.has(gene.uid)) {
          throw new Error("every gene needs a unique uid");
        }
        genes.add(gene.uid);
      }
    }
  }
  for (const group of data.groups) {
    if (!group || !Array.isArray(group.genes) || group.genes.some(uid => !genes.has(uid))) {
      throw new Error("every group must contain known gene IDs");
    }
  }
  for (const link of data.links) {
    if (!link?.query || !link?.target || !genes.has(link.query.uid) || !genes.has(link.target.uid)) {
      throw new Error("every link must refer to known genes");
    }
  }
}

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

function mountImportedPlot(data, { config = null, appearance = null, state = null } = {}) {
  validatePlotData(data);
  const chartConfig = config || appearance || chart.config();
  if (!chartConfig || typeof chartConfig !== "object" || Array.isArray(chartConfig)) {
    throw new Error("appearance settings must be an object");
  }
  replaceChart(chartConfig);
  plot.datum(data).call(chart);
  if (state) chart.state(state);
  mountSelectionToolbar();
  plotEditor.chart = chart;
  latestSimilarity = null;
  matrixButton.disabled = true;
  svgButton.disabled = false;
  dataButton.disabled = false;
  projectButton.disabled = false;
  editPlotButton.disabled = false;
  closePlotEditor();
  closeRunPanel();
}

async function readJsonFile(input) {
  const file = input.files[0];
  input.value = "";
  if (!file) return null;
  try {
    return JSON.parse(await file.text());
  } catch {
    throw new Error("the selected file is not valid JSON");
  }
}

projectButton.addEventListener("click", () => {
  const project = chart.project();
  const url = URL.createObjectURL(new Blob([`${JSON.stringify(project, null, 2)}\n`], {
    type: "application/json;charset=utf-8",
  }));
  const download = document.createElement("a");
  download.href = url;
  download.download = "clinker-project.json";
  document.body.append(download);
  download.click();
  download.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
});

loadDataButton.addEventListener("click", () => dataFileInput.click());
dataFileInput.addEventListener("change", async () => {
  try {
    const data = await readJsonFile(dataFileInput);
    if (!data) return;
    mountImportedPlot(data);
    setPlotStatus("Loaded plot data.");
  } catch (error) {
    setContextStatus(`Could not load plot data: ${error.message || String(error)}`, { tone: "error" });
  }
});

loadProjectButton.addEventListener("click", () => projectFileInput.click());
projectFileInput.addEventListener("change", async () => {
  try {
    const project = await readJsonFile(projectFileInput);
    if (!project) return;
    if (project.format !== "clinker-project" || project.version !== 1) {
      throw new Error("unsupported project format");
    }
    const projectConfig = project.config || project.appearance;
    if (!projectConfig) throw new Error("project is missing appearance settings");
    mountImportedPlot(project.data, {
      config: projectConfig,
      state: project.state,
    });
    setPlotStatus("Loaded clinker project.");
  } catch (error) {
    setContextStatus(`Could not load project: ${error.message || String(error)}`, { tone: "error" });
  }
});

function postProcess(layout, links, useFileOrder, onProgress) {
  const worker = new Worker(WORKER_URL, { type: "module" });
  return new Promise((resolve, reject) => {
    const finish = (callback, value) => {
      worker.terminate();
      callback(value);
    };
    worker.onmessage = ({ data }) => {
      if (data.type === "error") finish(reject, new Error(data.message));
      else if (data.type === "post-process-progress") onProgress(data);
      else if (data.type === "post-process") finish(resolve, data.result);
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage({
      type: "post-process",
      layout,
      useFileOrder,
      reportProgress: true,
      links: links.map(link => ({
        query: link.query,
        target: link.target,
        identity: link.identity,
        similarity: link.similarity,
      })),
    });
  });
}

function analyseTilesInWorkerPool(proteins, clusterCount, identity, candidatePairs, onProgress) {
  const clusters = proteinsByCluster(proteins, clusterCount);
  const totalTiles = candidatePairs
    ? Math.ceil(candidatePairs.length / PAIRS_PER_TILE)
    : tileCount(clusters);
  if (totalTiles === 0) return Promise.resolve({ links: [], alignmentCount: 0 });
  const tiles = candidatePairs ? candidateTiles(proteins, candidatePairs) : alignmentTiles(clusters);
  const totalPairs = candidatePairs ? candidatePairs.length : pairCount(clusters);
  const workers = Array.from(
    { length: workerCount(totalTiles) },
    () => new Worker(WORKER_URL, { type: "module" }),
  );
  const linksByTask = new Array(totalTiles);
  const alignmentCountsByTask = new Array(totalTiles).fill(0);
  const processedPairsByTask = new Array(totalTiles).fill(0);
  let nextTask = 0;
  let completed = 0;
  let settled = false;

  return new Promise((resolve, reject) => {
    const finish = error => {
      if (settled) return;
      settled = true;
      workers.forEach(worker => worker.terminate());
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

    const dispatch = worker => {
      if (nextTask === totalTiles) return;
      const taskIndex = nextTask++;
      const tile = tiles.next().value;
      processedPairsByTask[taskIndex] = 0;
      worker.postMessage({ ...tile, identity, taskIndex });
    };

    workers.forEach(worker => {
      worker.onmessage = ({ data }) => {
        if (settled) return;
        if (data.type === "error") {
          finish(new Error(data.message));
          return;
        }

        if (data.type === "tile-progress") {
          processedPairsByTask[data.taskIndex] = data.processedPairs;
          onProgress(processedPairsByTask.reduce((sum, count) => sum + count, 0), totalPairs);
          return;
        }

        linksByTask[data.taskIndex] = data.links.map(plotLink);
        alignmentCountsByTask[data.taskIndex] = data.alignmentCount || 0;
        processedPairsByTask[data.taskIndex] = Math.max(
          processedPairsByTask[data.taskIndex],
          data.pairCount || processedPairsByTask[data.taskIndex],
        );
        completed += 1;
        onProgress(processedPairsByTask.reduce((sum, count) => sum + count, 0), totalPairs);
        if (completed === totalTiles) finish();
        else dispatch(worker);
      };
      worker.onerror = event => finish(new Error(event.message || "Worker failed"));
      dispatch(worker);
    });
  });
}

async function analyseFiles(files) {
  const identity = Number(identityInput.value);
  if (!Number.isFinite(identity) || identity < 0 || identity > 1) {
    setRunStatus("Identity cutoff must be a number between 0 and 1.", { tone: "error" });
    return;
  }
  const kmerSize = Number(kmerSizeInput.value);
  const minSharedKmers = Number(minSharedKmersInput.value);
  if (prefilterInput.checked && (!Number.isInteger(kmerSize) || kmerSize < 1 || kmerSize > 8)) {
    setRunStatus("k-mer size must be a whole number between 1 and 8.", { tone: "error" });
    return;
  }
  if (prefilterInput.checked && (!Number.isInteger(minSharedKmers) || minSharedKmers < 1)) {
    setRunStatus("Minimum shared k-mers must be a positive whole number.", { tone: "error" });
    return;
  }

  analysisInProgress = true;
  setAnalysisInputsDisabled(true);
  matrixButton.disabled = true;
  svgButton.disabled = true;
  dataButton.disabled = true;
  projectButton.disabled = true;
  editPlotButton.disabled = true;
  closePlotEditor();
  let activeStep = "parse";
  try {
    const prefilter = {
      enabled: prefilterInput.checked,
      kmerSize,
      minSharedKmers,
      identityCutoff: identity,
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
        activeStep = "prefilter";
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
    });
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
    activeStep = "align";
    setRunStep("align", "active", "Aligning genes… 0%", { progress: 0 });
    const alignmentResult = await analyseTilesInWorkerPool(
      parsed.proteins,
      parsed.layout.clusters.length,
      identity,
      parsed.candidatePairs,
      (completedPairs, total) => {
        const progress = progressPercent(completedPairs, total);
        setRunStep("align", "active", `Aligning genes… ${progress}%`, { progress });
      },
    );
    const { links, alignmentCount } = alignmentResult;
    setRunStep(
      "align",
      "complete",
      `Ran ${alignmentCount.toLocaleString()} alignments · retained ${links.length.toLocaleString()} links`,
    );
    const useFileOrder = clusterOrderInput.value === "file";
    activeStep = "order";
    setRunStep("order", "active", "Scoring cluster pairs…", { progress: 0 });
    const result = await postProcess(
      parsed.layout,
      links,
      useFileOrder,
      ({ stage, completed, total }) => {
        if (stage === "similarity") {
          activeStep = "order";
          setRunStep(
            "order",
            "active",
            `Scored ${completed.toLocaleString()} of ${total.toLocaleString()} cluster pairs…`,
            { progress: progressPercent(completed, total) },
          );
        } else if (stage === "ordering") {
          activeStep = "order";
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
          activeStep = "groups";
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
          activeStep = "layout";
          setRunStep(
            "layout",
            "active",
            `Arranged ${completed.toLocaleString()} of ${total.toLocaleString()} clusters…`,
            { progress: progressPercent(completed, total) },
          );
        }
      },
    );
    const plotData = result.plotData;
    latestSimilarity = {
      clusterNames: result.clusterNames,
      clusterOrder: result.clusterOrder,
      similarityMatrix: result.similarityMatrix,
    };
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
    activeStep = null;

    replaceChart();
    plot.datum(plotData).call(chart);
    mountSelectionToolbar();
    plotEditor.chart = chart;
    svgButton.disabled = false;
    dataButton.disabled = false;
    projectButton.disabled = false;
    editPlotButton.disabled = false;
    analysisInProgress = false;
    closeRunPanel();
    setPlotStatus(summary);
  } catch (error) {
    if (activeStep) {
      const step = runSteps.get(activeStep);
      setRunStep(activeStep, "error", step.result.textContent);
    }
    setRunStatus(`Analysis failed: ${error.message || String(error)}`, {
      tone: "error",
      preserveSteps: true,
    });
  } finally {
    analysisInProgress = false;
    setAnalysisInputsDisabled(false);
    if (chart.data()) {
      matrixButton.disabled = !latestSimilarity;
      svgButton.disabled = false;
      dataButton.disabled = false;
      projectButton.disabled = false;
      editPlotButton.disabled = false;
    }
  }
}

analyseButton.addEventListener("click", async () => {
  if (filesInput.files.length === 0) {
    setRunStatus("Choose at least one GenBank file.", { tone: "error" });
    return;
  }
  try {
    setRunStatus("Reading selected files…");
    const files = await Promise.all([...filesInput.files].map(async file => ({
      name: file.name,
      bytes: new Uint8Array(await file.arrayBuffer()),
    })));
    await analyseFiles(files);
  } catch (error) {
    setRunStatus(`Analysis failed: ${error.message || String(error)}`, { tone: "error" });
  }
});

exampleButton.addEventListener("click", async () => {
  analyseButton.disabled = true;
  exampleButton.disabled = true;
  try {
    setRunStatus("Loading bundled example files…");
    const files = await Promise.all(EXAMPLE_FILES.map(async name => {
      const response = await fetch(`examples/${encodeURIComponent(name)}`);
      if (!response.ok) throw new Error(`Could not load ${name} (${response.status})`);
      return { name, bytes: new Uint8Array(await response.arrayBuffer()) };
    }));
    await analyseFiles(files);
  } catch (error) {
    setRunStatus(`Example analysis failed: ${error.message || String(error)}`, { tone: "error" });
  } finally {
    analyseButton.disabled = false;
    exampleButton.disabled = false;
  }
});
