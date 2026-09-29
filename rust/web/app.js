const filesInput = document.querySelector("#files");
const identityInput = document.querySelector("#identity");
const analyseButton = document.querySelector("#analyse");
const status = document.querySelector("#status");
const plot = d3.select("#plot");
let chart;

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

function clusterPairs(files) {
  return files.flatMap((_, queryIndex) =>
    files.slice(queryIndex + 1).map((_, offset) => [queryIndex, queryIndex + offset + 1]));
}

function parseClusterMetadata(files) {
  const worker = new Worker("worker.js", { type: "module" });
  const taskFiles = copyFiles(files, files.map((_, index) => index));

  return new Promise((resolve, reject) => {
    const finish = (callback, value) => {
      worker.terminate();
      callback(value);
    };
    worker.onmessage = ({ data }) => {
      if (data.type === "error") finish(reject, new Error(data.message));
      else finish(resolve, data.plotData);
    };
    worker.onerror = event => finish(reject, new Error(event.message || "Worker failed"));
    worker.postMessage(
      { type: "parse", files: taskFiles },
      taskFiles.map(file => file.bytes.buffer),
    );
  });
}

function plotLink(pairLink, [queryCluster, targetCluster]) {
  return {
    query: { uid: `gene-${queryCluster}-${pairLink.query.locus}-${pairLink.query.gene}` },
    target: { uid: `gene-${targetCluster}-${pairLink.target.locus}-${pairLink.target.gene}` },
    identity: pairLink.identity,
    similarity: pairLink.similarity,
  };
}

function analysePairsInWorkerPool(files, identity, onProgress) {
  const pairIndexes = clusterPairs(files);
  if (pairIndexes.length === 0) return Promise.resolve([]);
  const workers = Array.from(
    { length: workerCount(pairIndexes.length) },
    () => new Worker("worker.js", { type: "module" }),
  );
  const linksByTask = new Array(pairIndexes.length);
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
      if (nextTask === pairIndexes.length) return;
      const taskIndex = nextTask++;
      const fileIndexes = pairIndexes[taskIndex];
      const taskFiles = copyFiles(files, fileIndexes);
      worker.postMessage(
        { type: "analyse-pair", files: taskFiles, identity, fileIndexes, taskIndex },
        taskFiles.map(file => file.bytes.buffer),
      );
    };

    workers.forEach(worker => {
      worker.onmessage = ({ data }) => {
        if (settled) return;
        if (data.type === "error") {
          finish(new Error(data.message));
          return;
        }

        linksByTask[data.taskIndex] = data.links.map(link => plotLink(link, data.fileIndexes));
        completed += 1;
        onProgress(completed, pairIndexes.length);
        if (completed === pairIndexes.length) finish();
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
    const files = await Promise.all([...filesInput.files].map(async file => ({
      name: file.name,
      bytes: new Uint8Array(await file.arrayBuffer()),
    })));
    status.textContent = "Parsing GenBank files locally…";
    const metadata = await parseClusterMetadata(files);
    const links = await analysePairsInWorkerPool(files, identity, (completed, total) => {
      status.textContent = `Analysing locally… ${completed}/${total} cluster pairs`;
    });
    const plotData = { clusters: metadata.clusters, links, groups: [] };
    status.textContent = `${plotData.clusters.length} clusters; ${plotData.links.length} retained links.`;
    plot.selectAll("*").remove();
    chart = ClusterMap.ClusterMap();
    plot.datum(plotData).call(chart);
  } catch (error) {
    status.textContent = `Analysis failed: ${error.message || String(error)}`;
  } finally {
    analyseButton.disabled = false;
  }
});
