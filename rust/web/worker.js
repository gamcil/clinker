import init, { analyse_pair, parse_files } from "./pkg/clinker_wasm.js";

const wasm = init();

self.onmessage = async ({ data }) => {
  try {
    await wasm;
    if (data.type === "parse") {
      const plotData = parse_files(data.files);
      self.postMessage({ type: "parsed", plotData });
    } else if (data.type === "analyse-pair") {
      const links = analyse_pair(data.files, data.identity);
      self.postMessage({
        type: "pair-result",
        links,
        fileIndexes: data.fileIndexes,
        taskIndex: data.taskIndex,
      });
    }
  } catch (error) {
    self.postMessage({ type: "error", message: String(error) });
  }
};
