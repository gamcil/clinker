import init, {
  analyse_pairs,
  analyse_tile,
  BrowserFileParser,
  post_process,
} from "./pkg/clinker_wasm.js?v=20261006-wasm-init-v11";

const WASM_VERSION = "20261006-wasm-init-v11";
const wasmUrl = new URL("./pkg/clinker_wasm_bg.wasm", import.meta.url);
wasmUrl.searchParams.set("v", WASM_VERSION);
const wasm = init({ module_or_path: wasmUrl });
let fileParser = null;

function releaseFileParser() {
  if (!fileParser) return;
  fileParser.free();
  fileParser = null;
}

self.onmessage = async ({ data }) => {
  try {
    await wasm;
    if (data.type === "parse-start") {
      releaseFileParser();
      fileParser = new BrowserFileParser(data.prefilter);
      self.postMessage({ type: "parse-ready" });
    } else if (data.type === "parse-file") {
      if (!fileParser) throw new Error("file parser has not been started");
      fileParser.parse_file(data.file);
      self.postMessage({ type: "file-parsed", fileIndex: data.fileIndex });
    } else if (data.type === "parse-finish") {
      if (!fileParser) throw new Error("file parser has not been started");
      const filterProgress = progress => self.postMessage({
        type: "filter-progress",
        ...progress,
      });
      const parsed = fileParser.finish(filterProgress);
      releaseFileParser();
      self.postMessage({ type: "parsed", parsed });
    } else if (data.type === "parse") {
      // Compatibility for an already-open page running the previous app.js.
      releaseFileParser();
      fileParser = new BrowserFileParser(data.prefilter);
      for (let index = 0; index < data.files.length; index += 1) {
        fileParser.parse_file(data.files[index]);
        self.postMessage({
          type: "parse-progress",
          processedFiles: index + 1,
          totalFiles: data.files.length,
        });
        await new Promise(resolve => setTimeout(resolve, 0));
      }
      // The compatibility caller does not understand the new filtering
      // progress message and would mistake it for the final parse result.
      const parsed = fileParser.finish(() => {});
      releaseFileParser();
      self.postMessage({ type: "parsed", parsed });
    } else if (data.type === "analyse-tile") {
      const progress = processedPairs => self.postMessage({
        type: "tile-progress",
        taskIndex: data.taskIndex,
        processedPairs,
      });
      const result = analyse_tile(data.query, data.target, data.identity, progress);
      self.postMessage({
        type: "tile-result",
        links: result.links,
        alignmentCount: result.alignmentCount,
        taskIndex: data.taskIndex,
        pairCount: data.pairCount,
      });
    } else if (data.type === "analyse-pairs") {
      const progress = processedPairs => self.postMessage({
        type: "tile-progress",
        taskIndex: data.taskIndex,
        processedPairs,
      });
      const result = analyse_pairs(data.proteins, data.pairs, data.identity, progress);
      self.postMessage({
        type: "tile-result",
        links: result.links,
        alignmentCount: result.alignmentCount,
        taskIndex: data.taskIndex,
        pairCount: data.pairCount,
      });
    } else if (data.type === "post-process") {
      const progress = data.reportProgress
        ? update => self.postMessage({ type: "post-process-progress", ...update })
        : () => {};
      const result = post_process(data.layout, data.links, data.useFileOrder, progress);
      self.postMessage({ type: "post-process", result });
    }
  } catch (error) {
    releaseFileParser();
    self.postMessage({ type: "error", message: String(error) });
  }
};
