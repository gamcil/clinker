import init, { analyse_pairs, analyse_tile, parse_files, post_process } from "./pkg/clinker_wasm.js";

const wasm = init();

self.onmessage = async ({ data }) => {
  try {
    await wasm;
    if (data.type === "parse") {
      const parsed = parse_files(data.files, data.prefilter);
      self.postMessage({ type: "parsed", parsed });
    } else if (data.type === "analyse-tile") {
      const links = analyse_tile(data.query, data.target, data.identity);
      self.postMessage({
        type: "tile-result",
        links,
        taskIndex: data.taskIndex,
      });
    } else if (data.type === "analyse-pairs") {
      const links = analyse_pairs(data.proteins, data.pairs, data.identity);
      self.postMessage({
        type: "tile-result",
        links,
        taskIndex: data.taskIndex,
      });
    } else if (data.type === "post-process") {
      const result = post_process(data.layout, data.links);
      self.postMessage({ type: "post-process", result });
    }
  } catch (error) {
    self.postMessage({ type: "error", message: String(error) });
  }
};
