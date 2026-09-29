import init, { analyse } from "./pkg/clinker_wasm.js";

const wasm = init();

self.onmessage = async ({ data }) => {
  if (data.type !== "analyse") return;
  try {
    await wasm;
    const plotData = analyse(data.files, data.identity);
    self.postMessage({ type: "result", plotData });
  } catch (error) {
    self.postMessage({ type: "error", message: String(error) });
  }
};
