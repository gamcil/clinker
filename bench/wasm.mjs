import { Workflow } from '/rust/web/workflow.mjs';
import { analyseTilesInWorkerPool } from '/rust/web/alignment.mjs';

const workerUrl = '/rust/web/worker.js';
const compiled = await WebAssembly.compileStreaming(fetch('/rust/web/pkg/clinker_wasm_bg.wasm'));

function request(worker, message, expected, transfer = []) {
  return new Promise((resolve, reject) => {
    worker.onmessage = ({ data }) => {
      if (data.type === 'error') reject(new Error(data.message));
      else if (data.type === expected) resolve(data);
    };
    worker.onerror = event => reject(new Error(event.message));
    worker.onmessageerror = () => reject(new Error('Invalid worker response'));
    worker.postMessage(message, transfer);
  });
}

window.runBenchmark = async options => {
  // Fetch/file staging and module compilation are outside processing timers.
  const files = [];
  for (let index = 0; index < options.files.length; index++) {
    const response = await fetch(`/input/${index}`);
    if (!response.ok) throw new Error(`Input ${index}: HTTP ${response.status}`);
    files.push({ name: options.files[index].split('/').pop(), bytes: new Uint8Array(await response.arrayBuffer()) });
  }
  const workflow = new Workflow();
  try {
    const parser = workflow.worker(workerUrl);
    await request(parser, { type: 'parse-start', wasmModule: compiled, prefilter: {
      enabled: options.prefilter, kmerSize: options.kmer_size, minSharedKmers: options.min_shared_kmers,
      identityCutoff: options.identity,
    } }, 'parse-ready');
    const start = performance.now();
    for (let index = 0; index < files.length; index++) {
      await request(parser, { type: 'parse-file', file: files[index], fileIndex: index }, 'file-parsed', [files[index].bytes.buffer]);
    }
    const parsedAt = performance.now();
    const { parsed } = await request(parser, { type: 'parse-finish' }, 'parsed');
    const preparedAt = performance.now();
    workflow.releaseWorker(parser);
    if (options.inspect) {
      const proteins = parsed.layout.clusters.map(() => []);
      for (const protein of parsed.proteins) proteins[protein.cluster].push(protein.translation);
      return { proteins };
    }
    const alignAt = performance.now();
    const result = await analyseTilesInWorkerPool(
      parsed, options.identity, compiled, () => {}, workflow, options.workers,
    );
    const end = performance.now();
    return {
      parse_ms: parsedAt - start,
      prepare_ms: options.prefilter ? 0 : preparedAt - parsedAt,
      prefilter_ms: options.prefilter ? preparedAt - parsedAt : null,
      alignment_ms: end - alignAt, total_ms: end - start,
      clusters: parsed.layout.clusters.length, proteins: parsed.proteins.length,
      amino_acids: parsed.proteins.reduce((n, p) => n + p.translation.length, 0),
      possible_pairs: parsed.alignmentCounts.possiblePairs,
      eligible_pairs: null,
      candidate_pairs: parsed.alignmentCounts.plannedPairs,
      aligned_pairs: result.alignmentCount, retained_links: result.links.length,
    };
  } finally { workflow.dispose(); }
};
