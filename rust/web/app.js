const filesInput = document.querySelector("#files");
const identityInput = document.querySelector("#identity");
const prefilterInput = document.querySelector("#prefilter");
const analyseButton = document.querySelector("#analyse");
const exampleButton = document.querySelector("#load-example");
const matrixButton = document.querySelector("#download-matrix");
const svgButton = document.querySelector("#download-svg");
const status = document.querySelector("#status");
const plot = d3.select("#plot");
const GENES_PER_TILE_SIDE = 40;
const PAIRS_PER_TILE = 1600;

function createChart() {
  return ClusterMap.ClusterMap().config({
    link: { bestOnly: true },
    plot: { renderer: "webgpu" },
    legend: { columns: 4, position: "bottom" },
  });
}

// A chart keeps its camera and editable locus state. Each completed analysis
// is a new dataset, so give it a fresh instance rather than carrying state
// between positional IDs such as `cluster-0` and `locus-0-0`.
let chart = createChart();
let latestSimilarity = null;
const EXAMPLE_FILES = [
  "P. vexata CBS 129021.gbk",
  "A. versicolor CBS 583.65.gbk",
  "A. mulundensis DSM 5745.gbk",
  "A. burnettii MST-FP2249.gbk",
  "A. alliaceus CBS 536.65.gbk",
];

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

function parseInputFiles(files, prefilter) {
  const worker = new Worker("worker.js", { type: "module" });
  const taskFiles = copyFiles(files, files.map((_, index) => index));

  return new Promise((resolve, reject) => {
    const finish = (callback, value) => {
      worker.terminate();
      callback(value);
    };
    worker.onmessage = ({ data }) => {
      if (data.type === "error") finish(reject, new Error(data.message));
      else finish(resolve, data.parsed);
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage(
      { type: "parse", files: taskFiles, prefilter },
      taskFiles.map(file => file.bytes.buffer),
    );
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

function bestLinksForDisplay(plotData) {
  // The pinned clustermap refactor applies `link.bestOnly` in SVG but not yet
  // in its WebGPU renderer. Keep that renderer workaround local: Rust retains
  // every link for grouping and output, while this copy contains only the
  // highest-identity overlapping link for each cluster pair.
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

  return [...selectedByPair.values()].flat();
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
    status.textContent = `SVG export failed: ${error.message || String(error)}`;
  }
});

function postProcess(layout, links) {
  const worker = new Worker("worker.js", { type: "module" });
  return new Promise((resolve, reject) => {
    const finish = (callback, value) => {
      worker.terminate();
      callback(value);
    };
    worker.onmessage = ({ data }) => {
      if (data.type === "error") finish(reject, new Error(data.message));
      else finish(resolve, data.result);
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage({
      type: "post-process",
      layout,
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
  if (totalTiles === 0) return Promise.resolve([]);
  const tiles = candidatePairs ? candidateTiles(proteins, candidatePairs) : alignmentTiles(clusters);
  const totalPairs = candidatePairs ? candidatePairs.length : pairCount(clusters);
  const workers = Array.from(
    { length: workerCount(totalTiles) },
    () => new Worker("worker.js", { type: "module" }),
  );
  const linksByTask = new Array(totalTiles);
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
        resolve(links);
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
    status.textContent = "Identity cutoff must be a number between 0 and 1.";
    return;
  }

  analyseButton.disabled = true;
  exampleButton.disabled = true;
  matrixButton.disabled = true;
  svgButton.disabled = true;
  latestSimilarity = null;
  try {
    status.textContent = "Parsing GenBank files…";
    const prefilter = {
      enabled: prefilterInput.checked,
      kmerSize: 3,
      minSharedKmers: 3,
      identityCutoff: identity,
    };
    const parsed = await parseInputFiles(files, prefilter);
    status.textContent = "Aligning genes… 0%";
    const links = await analyseTilesInWorkerPool(
      parsed.proteins,
      parsed.layout.clusters.length,
      identity,
      parsed.candidatePairs,
      (completedPairs, total) => {
        status.textContent = `Aligning genes… ${progressPercent(completedPairs, total)}%`;
      },
    );
    status.textContent = "Building homology groups and ordering clusters…";
    const result = await postProcess(parsed.layout, links);
    const plotData = result.plotData;
    latestSimilarity = {
      clusterNames: result.clusterNames,
      clusterOrder: result.clusterOrder,
      similarityMatrix: result.similarityMatrix,
    };
    matrixButton.disabled = false;
    const displayData = {
      ...plotData,
      links: bestLinksForDisplay(plotData),
    };
    status.textContent = `${plotData.clusters.length} clusters; ${plotData.links.length} retained links; ${displayData.links.length} displayed best links.`;

    chart.destroy();
    chart = createChart();
    plot.datum(displayData).call(chart);
    svgButton.disabled = false;
  } catch (error) {
    status.textContent = `Analysis failed: ${error.message || String(error)}`;
  } finally {
    analyseButton.disabled = false;
    exampleButton.disabled = false;
  }
}

analyseButton.addEventListener("click", async () => {
  if (filesInput.files.length === 0) {
    status.textContent = "Choose at least one GenBank file.";
    return;
  }
  try {
    status.textContent = "Reading files…";
    const files = await Promise.all([...filesInput.files].map(async file => ({
      name: file.name,
      bytes: new Uint8Array(await file.arrayBuffer()),
    })));
    await analyseFiles(files);
  } catch (error) {
    status.textContent = `Analysis failed: ${error.message || String(error)}`;
  }
});

exampleButton.addEventListener("click", async () => {
  analyseButton.disabled = true;
  exampleButton.disabled = true;
  try {
    status.textContent = "Loading bundled example files…";
    const files = await Promise.all(EXAMPLE_FILES.map(async name => {
      const response = await fetch(`examples/${encodeURIComponent(name)}`);
      if (!response.ok) throw new Error(`Could not load ${name} (${response.status})`);
      return { name, bytes: new Uint8Array(await response.arrayBuffer()) };
    }));
    await analyseFiles(files);
  } catch (error) {
    status.textContent = `Example analysis failed: ${error.message || String(error)}`;
  } finally {
    analyseButton.disabled = false;
    exampleButton.disabled = false;
  }
});
