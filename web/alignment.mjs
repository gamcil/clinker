const WORKER_URL = new URL("./worker.js?v=20261020-shared-preparation-v1", import.meta.url);

function workerCount(taskCount, requestedWorkers) {
  if (requestedWorkers !== undefined) {
    if (!Number.isInteger(requestedWorkers) || requestedWorkers < 1) throw new Error('Worker count must be a positive integer');
    return Math.min(taskCount, requestedWorkers);
  }
  // Each worker owns a separate Wasm instance and its alignment buffers. A
  // modest cap avoids making many large DP matrices resident at once.
  return Math.min(taskCount, Math.max(1, Math.min(navigator.hardwareConcurrency || 2, 4)));
}

// Rust plans the batches; JavaScript only packs their proteins for transfer.
function* alignmentTiles(proteins, candidatePairs, batches) {
  for (const batch of batches) {
    if (batch.kind === "cartesian") {
      const query = proteins.slice(batch.query.start, batch.query.end);
      const target = proteins.slice(batch.target.start, batch.target.end);
      yield {
        proteins: [...query, ...target], pairs: [],
        batch: { kind: "cartesian", query: { start: 0, end: query.length }, target: { start: query.length, end: query.length + target.length } },
      };
    } else {
      const proteinIndexes = new Map();
      const tileProteins = [];
      const localIndex = index => {
        if (!proteinIndexes.has(index)) {
          proteinIndexes.set(index, tileProteins.length);
          tileProteins.push(proteins[index]);
        }
        return proteinIndexes.get(index);
      };
      const pairs = candidatePairs.slice(batch.pairs.start, batch.pairs.end).map(pair => ({
        queryIndex: localIndex(pair.queryIndex),
        targetIndex: localIndex(pair.targetIndex),
      }));
      yield { proteins: tileProteins, pairs, batch: { kind: "selected", pairs: { start: 0, end: pairs.length } } };
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

function analyseTilesInWorkerPool(
  prepared,
  identity,
  wasmModule,
  onProgress,
  workflow,
  requestedWorkers,
) {
  const { proteins, candidatePairs, alignmentBatches: batches, alignmentCounts: counts } = prepared;
  const totalTiles = counts.batchCount;
  if (totalTiles === 0) return Promise.resolve({ links: [], alignmentCount: 0 });
  const tiles = alignmentTiles(proteins, candidatePairs, batches);
  const totalPairs = counts.plannedPairs;
  const workers = Array.from(
    { length: workerCount(totalTiles, requestedWorkers) },
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
      if (settled || nextTask === totalTiles) return;
      const taskIndex = nextTask++;
      try {
        const tile = tiles.next().value;
        processedPairsByTask[taskIndex] = 0;
        worker.postMessage({ type: "analyse-batch", ...tile, identity, taskIndex, wasmModule });
      } catch (error) { finish(error); }
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
      worker.onmessageerror = () => finish(new Error("Could not read the analysis worker response."));
      dispatch(worker);
    });
  });
}


export { analyseTilesInWorkerPool };
