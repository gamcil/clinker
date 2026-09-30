const filesInput = document.querySelector("#files");
const identityInput = document.querySelector("#identity");
const analyseButton = document.querySelector("#analyse");
const status = document.querySelector("#status");
const plot = d3.select("#plot");
const GENES_PER_TILE_SIDE = 20;
// Keep one chart instance, as the original clustermap integration does. The
// library retains its renderer state on this object between redraws.
const chart = ClusterMap.ClusterMap();

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

function parseInputFiles(files) {
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
      { type: "parse", files: taskFiles },
      taskFiles.map(file => file.bytes.buffer),
    );
  });
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

function* alignmentTiles(clusters) {
  for (let queryCluster = 0; queryCluster < clusters.length; queryCluster += 1) {
    for (let targetCluster = queryCluster + 1; targetCluster < clusters.length; targetCluster += 1) {
      const queryProteins = clusters[queryCluster];
      const targetProteins = clusters[targetCluster];
      for (let queryStart = 0; queryStart < queryProteins.length; queryStart += GENES_PER_TILE_SIDE) {
        for (let targetStart = 0; targetStart < targetProteins.length; targetStart += GENES_PER_TILE_SIDE) {
          yield {
            query: queryProteins.slice(queryStart, queryStart + GENES_PER_TILE_SIDE),
            target: targetProteins.slice(targetStart, targetStart + GENES_PER_TILE_SIDE),
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

function analyseTilesInWorkerPool(proteins, clusterCount, identity, onProgress) {
  const clusters = proteinsByCluster(proteins, clusterCount);
  const totalTiles = tileCount(clusters);
  if (totalTiles === 0) return Promise.resolve([]);
  const tiles = alignmentTiles(clusters);
  const workers = Array.from(
    { length: workerCount(totalTiles) },
    () => new Worker("worker.js", { type: "module" }),
  );
  const linksByTask = new Array(totalTiles);
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
      worker.postMessage({ type: "analyse-tile", ...tile, identity, taskIndex });
    };

    workers.forEach(worker => {
      worker.onmessage = ({ data }) => {
        if (settled) return;
        if (data.type === "error") {
          finish(new Error(data.message));
          return;
        }

        linksByTask[data.taskIndex] = data.links.map(plotLink);
        completed += 1;
        onProgress(completed, totalTiles);
        if (completed === totalTiles) finish();
        else dispatch(worker);
      };
      worker.onerror = event => finish(new Error(event.message || "Worker failed"));
      dispatch(worker);
    });
  });
}

analyseButton.addEventListener("click", async () => {
  if (filesInput.files.length === 0) {
    status.textContent = "Choose at least one GenBank file.";
    return;
  }

  const identity = Number(identityInput.value);
  if (!Number.isFinite(identity) || identity < 0 || identity > 1) {
    status.textContent = "Identity cutoff must be a number between 0 and 1.";
    return;
  }

  analyseButton.disabled = true;
  try {
    status.textContent = "Reading files…";
    let files = await Promise.all([...filesInput.files].map(async file => ({
      name: file.name,
      bytes: new Uint8Array(await file.arrayBuffer()),
    })));
    status.textContent = "Parsing GenBank files locally…";
    const parsed = await parseInputFiles(files);
    files = null;
    const totalTiles = tileCount(proteinsByCluster(parsed.proteins, parsed.layout.clusters.length));
    status.textContent = `Analysing locally… 0/${totalTiles} protein tiles`;
    const links = await analyseTilesInWorkerPool(
      parsed.proteins,
      parsed.layout.clusters.length,
      identity,
      (completed, total) => {
        status.textContent = `Analysing locally… ${completed}/${total} protein tiles`;
      },
    );
    status.textContent = "Building homology groups and ordering clusters…";
    const result = await postProcess(parsed.layout, links);
    const plotData = result.plotData;
    status.textContent = `${plotData.clusters.length} clusters; ${plotData.links.length} retained links.`;

    // clustermap animates parts of a render. Stop those transitions before
    // removing its SVG so an old render cannot try to update detached nodes.
    plot.selectAll("*").interrupt();
    plot.selectAll("*").remove();
    // ClusterMap expects its input through a one-element D3 data join (the
    // same pattern used by clinker/plot/clinker.js), rather than `datum`.
    plot.data([plotData]).call(chart);
  } catch (error) {
    status.textContent = `Analysis failed: ${error.message || String(error)}`;
  } finally {
    analyseButton.disabled = false;
  }
});
