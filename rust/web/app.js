const filesInput = document.querySelector("#files");
const identityInput = document.querySelector("#identity");
const analyseButton = document.querySelector("#analyse");
const status = document.querySelector("#status");
const plot = d3.select("#plot");
const worker = new Worker("worker.js", { type: "module" });

let chart;

worker.onmessage = ({ data }) => {
  analyseButton.disabled = false;
  if (data.type === "error") {
    status.textContent = `Analysis failed: ${data.message}`;
    return;
  }

  const { plotData } = data;
  status.textContent = `${plotData.clusters.length} clusters; ${plotData.links.length} retained links.`;
  plot.selectAll("*").remove();
  chart = ClusterMap.ClusterMap();
  plot.datum(plotData).call(chart);
};

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
  status.textContent = "Reading files…";
  const files = await Promise.all([...filesInput.files].map(async file => ({
    name: file.name,
    bytes: new Uint8Array(await file.arrayBuffer()),
  })));
  const transfer = files.map(file => file.bytes.buffer);
  status.textContent = "Analysing locally…";
  worker.postMessage({ type: "analyse", files, identity }, transfer);
});
