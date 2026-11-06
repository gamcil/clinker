import { Worker, isMainThread, parentPort, workerData } from "node:worker_threads";
import { readFile, readdir, writeFile } from "node:fs/promises";
import { fileURLToPath } from "node:url";
import path from "node:path";
import { performance } from "node:perf_hooks";

const packageUrl = new URL("../../web/.build/pkg/clinker_wasm.js", import.meta.url);
const wasmUrl = new URL("../../web/.build/pkg/clinker_wasm_bg.wasm", import.meta.url);

if (!isMainThread) {
  const { default: init, AlignmentWorker } = await import(packageUrl.href);
  await init({ module_or_path: await readFile(wasmUrl) });
  parentPort.postMessage({ ready: true });
  const aligner = new AlignmentWorker();
  parentPort.on("message", tile => {
    try {
      const result = aligner.analyse_batch(tile.proteins, tile.batch, tile.pairs, workerData.identity, () => {});
      parentPort.postMessage({ links: result.links.length });
    } catch (error) {
      parentPort.postMessage({ error: String(error) });
    }
  });
} else {
  const [dataset, identity, jobs, output] = process.argv.slice(2);
  const { default: init, BrowserFileParser } = await import(packageUrl.href);
  await init({ module_or_path: await readFile(wasmUrl) });

  const parser = new BrowserFileParser({
    enabled: false,
    noAlign: false,
    kmerSize: 5,
    minSharedKmers: 8,
    identityCutoff: Number(identity),
    genbank: {},
  });
  for (const name of (await readdir(dataset)).filter(name => name.endsWith(".gbk")).sort()) {
    const bytes = await readFile(path.join(dataset, name));
    parser.parse_file({ name, bytes: new Uint8Array(bytes) });
  }

  const workerCount = Math.max(1, Math.min(Number(jobs), 4));
  const workers = Array.from({ length: workerCount }, () => new Worker(fileURLToPath(import.meta.url), {
    workerData: { identity: Number(identity) },
  }));
  await Promise.all(workers.map(worker => new Promise((resolve, reject) => {
    worker.once("message", message => message.ready
      ? resolve()
      : reject(new Error("WASM alignment worker failed to initialize")));
    worker.once("error", reject);
  })));
  let nextBatch = 0;
  let completed = 0;
  let retainedLinks = 0;
  let settled = false;
  let dispatch;
  const started = performance.now();
  const prepared = parser.finish(() => {});
  const tiles = alignmentTiles(prepared);
  const finish = async error => {
    if (settled) return;
    settled = true;
    await Promise.all(workers.map(worker => worker.terminate()));
    if (error) throw error;
    const elapsedSeconds = (performance.now() - started) / 1000;
    await writeFile(output, `${JSON.stringify({ elapsed_seconds: elapsedSeconds, retained_links: retainedLinks })}\n`);
  };
  const done = new Promise((resolve, reject) => {
    dispatch = worker => {
      if (nextBatch === prepared.alignmentBatches.length) return;
      worker.postMessage(tiles.next().value);
      nextBatch += 1;
    };
    for (const worker of workers) {
      worker.on("message", async result => {
        if (result.error) {
          try { await finish(new Error(result.error)); reject(new Error(result.error)); } catch (error) { reject(error); }
          return;
        }
        retainedLinks += result.links;
        completed += 1;
        if (completed === prepared.alignmentBatches.length) {
          try { await finish(); resolve(); } catch (error) { reject(error); }
        } else {
          dispatch(worker);
        }
      });
      worker.on("error", reject);
      dispatch(worker);
    }
  });
  if (prepared.alignmentBatches.length === 0) {
    await finish();
    process.exit(0);
  }
  await done;
}

function* alignmentTiles(prepared) {
  const { proteins, candidatePairs, alignmentBatches } = prepared;
  for (const batch of alignmentBatches) {
    if (batch.kind === "cartesian") {
      const query = proteins.slice(batch.query.start, batch.query.end);
      const target = proteins.slice(batch.target.start, batch.target.end);
      yield {
        proteins: [...query, ...target],
        pairs: [],
        batch: {
          kind: "cartesian",
          query: { start: 0, end: query.length },
          target: { start: query.length, end: query.length + target.length },
        },
      };
    } else {
      const indexes = new Map();
      const tileProteins = [];
      const localIndex = index => {
        if (!indexes.has(index)) {
          indexes.set(index, tileProteins.length);
          tileProteins.push(proteins[index]);
        }
        return indexes.get(index);
      };
      const pairs = candidatePairs.slice(batch.pairs.start, batch.pairs.end).map(pair => ({
        queryIndex: localIndex(pair.queryIndex),
        targetIndex: localIndex(pair.targetIndex),
      }));
      yield {
        proteins: tileProteins,
        pairs,
        batch: { kind: "selected", pairs: { start: 0, end: pairs.length } },
      };
    }
  }
}
