import { mountRunInputs } from "./run-inputs.mjs";
import { parseInputFiles, postProcess, openAnalysisSession } from "./analysis.mjs";
import { analyseTilesInWorkerPool } from "./alignment.mjs";
import { parseProject } from "./viewer/project.mjs";
import { Workflow } from "./workflow.mjs";
import { mountTreeImport } from "./tree-import.mjs";
import { mountAnnotationImport } from "./annotation-import.mjs";
import { bindInputTree } from "./viewer/tree.mjs";
import { createViewer } from "./viewer/viewer.mjs";
import { readThemePreference, mountViewerControls } from "./viewer/controls.mjs";

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
const annotationImportButton = document.querySelector("#import-annotations");
const addDataMenu = document.querySelector("#add-data-menu");
const treeButton = document.querySelector("#manage-tree");
const plotOrder = document.querySelector("#plot-order");
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
const plotElement = document.querySelector("#plot");
const plotStage = document.querySelector(".plot-stage");
const viewer = createViewer({ plot: plotElement, stage: plotStage, editor: plotEditor, appearance: readThemePreference() });
let latestSimilarity = null;
let activeWorkflow = null;
let treeImport = null;
const runInputs = mountRunInputs(runPanel, {
  onChange: syncControls,
  onInputsChanged: () => treeImport?.inputsChanged(),
  onError: message => setRunStatus(message, { tone: "error" }),
});
const EXAMPLE_FILES = [
  "P. vexata CBS 129021.gbk",
  "A. versicolor CBS 583.65.gbk",
  "A. mulundensis DSM 5745.gbk",
  "A. burnettii MST-FP2249.gbk",
  "A. alliaceus CBS 536.65.gbk",
];

function syncControls() {
  const busy = activeWorkflow !== null;
  const hasPlot = Boolean(viewer.chart.data());
  runInputs.setDisabled(busy);
  analyseButton.disabled = busy || !runInputs.ready || !(treeImport?.runReady() ?? true);
  for (const button of [exampleButton, loadClinkerFileButton, openProjectStartButton, newRunButton]) button.disabled = busy;
  matrixButton.disabled = busy || !hasPlot || !latestSimilarity;
  svgButton.disabled = busy || !hasPlot;
  projectButton.disabled = busy || !hasPlot || !viewer.analysis;
  editPlotButton.disabled = busy || !hasPlot;
  treeButton.disabled = busy || !hasPlot;
  annotationImportButton.disabled = busy || !hasPlot;
  plotActions.hidden = !hasPlot;
  editPlotButton.hidden = !hasPlot;
  newRunButton.hidden = !hasPlot;
  newRunButton.setAttribute("aria-expanded", String(!runPanel.hidden));
  plotActions.inert = busy;
  plotOrder.inert = busy || !hasPlot;
  if (busy) { addDataMenu.open = false; plotOrder.open = false; }
  cancelRunButton.hidden = !busy && (!hasPlot || runPanel.hidden);
  cancelRunButton.textContent = busy ? "Cancel" : "Back to plot";
  document.querySelector(".run-card").classList.toggle("is-busy", busy);
}

async function runWorkflow(task, failureMessage) {
  if (activeWorkflow) return;
  treeImport.close();
  const workflow = new Workflow();
  const trigger = document.activeElement;
  openRunPanel();
  activeWorkflow = workflow;
  syncControls();
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
    syncControls();
    if (completed) {
      closeRunPanel();
      editPlotButton.focus();
    } else if (trigger instanceof HTMLElement && !trigger.disabled && trigger.getClientRects().length) {
      trigger.focus();
    }
  }
}

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

mountViewerControls({
  viewer, editor: plotEditor,
  getTrees: () => treeImport.projectTrees(),
  getSimilarity: () => latestSimilarity,
  onEditorOpen: closeRunPanel,
  onError: setContextStatus,
});

document.addEventListener("keydown", event => {
  if (event.key !== "Escape" || event.defaultPrevented) return;
  if (document.querySelector("dialog[open]")) return;
  if (!activeWorkflow && !runPanel.hidden && viewer.chart.data()) closeRunPanel();
});

function closeRunPanel() {
  if (activeWorkflow || !viewer.chart.data()) return;
  runPanel.hidden = true;
  syncControls();
}

function openRunPanel() {
  runPanel.hidden = false;
  setRunStatus("");
  viewer.closeEditor();
  syncControls();
}

newRunButton.addEventListener("click", () => {
  runInputs.clear();
  openRunPanel();
  runInputs.focus();
});
cancelRunButton.addEventListener("click", () => {
  if (activeWorkflow) activeWorkflow.cancel();
  else closeRunPanel();
});

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

function mountImportedPlot(saved) {
  treeImport.close();
  viewer.loadProject(saved);
  treeImport.bind();
  viewer.closeEditor();
}

function retainAnalysisResult(result) {
  latestSimilarity = result.similarityMatrix;
  viewer.setAnalysis(result.analysisSession);
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
    mountImportedPlot(saved);
    treeImport.restoreTrees(saved);
    retainAnalysisResult(analysisResult);
    setPlotStatus("Loaded clinker project with its analysis data.");
    return true;
  }, "Could not open project");
});

async function analyseFiles(files, workflow) {
  const runTree = treeImport.runTree();
  const { settings, useFileOrder, groupDefinitions } = await runInputs.readOptions(workflow);
  workflow.throwIfCancelled();
  const identity = settings.identityCutoff;
  showRunSteps({ includePrefilter: settings.enabled });
  setRunStep("parse", "active", `Parsing 0 of ${files.length} files…`, { progress: 0 });
  const parsed = await parseInputFiles(files, settings, {
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
  const resolvedRunTree = await treeImport.resolveRunTree(parsed.layout.clusters, runTree, settings.genbank.asSeparateClusters, workflow);
  workflow.throwIfCancelled();
  setRunStep(
    "parse",
    "complete",
    `Parsed ${files.length} ${files.length === 1 ? "file" : "files"} · ${parsed.proteins.length.toLocaleString()} genes`,
  );
  const { possiblePairs: possibleAlignments, plannedPairs: plannedAlignments } = parsed.alignmentCounts;
  if (settings.enabled) {
    setRunStep(
      "prefilter",
      "complete",
      prefilterResult(plannedAlignments, possibleAlignments),
    );
  }
  setRunStep("align", "active", "Aligning genes… 0%", { progress: 0 });
  const alignmentResult = settings.noAlign
    ? { links: [], alignmentCount: 0 }
    : await analyseTilesInWorkerPool(
      parsed,
      identity,
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
    settings.noAlign
      ? "Skipped protein alignment"
      : `Ran ${alignmentCount.toLocaleString()} alignments · retained ${links.length.toLocaleString()} links`,
  );
  setRunStep("order", "active", "Scoring cluster pairs…", { progress: 0 });
  const result = await postProcess(
    { layout: parsed.layout, proteins: parsed.proteins, links, groupDefinitions,
      settings, useFileOrder, wasmModule: parsed.wasmModule },
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
  if (resolvedRunTree) plotData.tree = bindInputTree(resolvedRunTree, plotData.clusters, settings.genbank.asSeparateClusters);
  const summary = `${plotData.clusters.length} clusters; ${plotData.links.length} retained links.`;
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

  treeImport.close();
  viewer.mount(plotData, { analysis: result.analysisSession });
  latestSimilarity = result.similarityMatrix;
  treeImport.bind();
  treeImport.setDendrogram(result.clusteringDendrogram);
  viewer.refreshEditor();
  setPlotStatus(summary);
  return true;
}

analyseButton.addEventListener("click", () => {
  const selected = runInputs.files;
  if (selected.length === 0) return;
  runWorkflow(async workflow => {
    setRunStatus("Preparing selected files…");
    const validation = await workflow.wait(runInputs.inspect());
    if (!validation.ready) throw new Error("Resolve the selected input-file problems before starting the analysis.");
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

treeImport = mountTreeImport({
  button: treeButton,
  getChart: () => viewer.chart,
  getRunInputs: () => runInputs.treeInputs(),
  onAppearance: () => {
    return viewer.openTreeAppearance();
  },
  onRunChange: syncControls,
  onChange: message => {
    viewer.refreshEditor();
    setPlotStatus(message);
  },
});

mountAnnotationImport({
  button: annotationImportButton,
  getChart: () => viewer.chart,
  onImported: message => setPlotStatus(message),
});

syncControls();
