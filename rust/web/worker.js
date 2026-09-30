import init, { analyse_tile, build_groups, parse_files } from "./pkg/clinker_wasm.js";

const wasm = init();

self.onmessage = async ({ data }) => {
  try {
    await wasm;
    if (data.type === "parse") {
      const parsed = parse_files(data.files);
      self.postMessage({ type: "parsed", parsed });
    } else if (data.type === "analyse-tile") {
      const links = analyse_tile(data.query, data.target, data.identity);
      self.postMessage({
        type: "tile-result",
        links,
        taskIndex: data.taskIndex,
      });
    } else if (data.type === "groups") {
      const groups = build_groups(data.links);
      self.postMessage({ type: "groups", groups });
    }
  } catch (error) {
    self.postMessage({ type: "error", message: String(error) });
  }
};
