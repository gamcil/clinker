export function displayedClusterOrder({ data, state }) {
  return state?.clusterOrder?.length ? state.clusterOrder : data.clusters.map(cluster => cluster.uid);
}

function csvCell(value) {
  const text = String(value);
  return /[",\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export function similarityCsv({ clusterNames, similarityMatrix }, clusterOrder) {
  const plotPosition = new Map(clusterOrder.map((cluster, index) => [cluster, index + 1]));
  if (clusterNames.some((_, index) => !plotPosition.has(`cluster-${index}`))) {
    throw new Error("Figure and analysis cluster IDs do not match.");
  }
  const rows = [[
    "query_plot_order",
    "query_cluster",
    "query_index",
    "target_plot_order",
    "target_cluster",
    "target_index",
    "query_to_target_coverage",
    "target_to_query_coverage",
    "containment_similarity",
    "distance",
  ]];
  for (let query = 0; query < clusterNames.length; query += 1) {
    for (let target = query + 1; target < clusterNames.length; target += 1) {
      const pair = similarityMatrix[query][target];
      rows.push([
        plotPosition.get(`cluster-${query}`),
        clusterNames[query],
        query,
        plotPosition.get(`cluster-${target}`),
        clusterNames[target],
        target,
        pair.queryCoverage,
        pair.targetCoverage,
        pair.similarity,
        1 - pair.similarity,
      ]);
    }
  }
  return `${rows.map(row => row.map(csvCell).join(",")).join("\n")}\n`;
}

function validatePlotData(data) {
  if (!data || typeof data !== "object" || Array.isArray(data)) {
    throw new Error("the file must contain a plot-data object");
  }
  if (!Array.isArray(data.clusters) || !Array.isArray(data.links) || !Array.isArray(data.groups)) {
    throw new Error("plot data must contain clusters, links, and groups arrays");
  }
  const genes = new Set();
  for (const cluster of data.clusters) {
    if (!cluster || !Array.isArray(cluster.loci)) throw new Error("every cluster needs a loci array");
    for (const locus of cluster.loci) {
      if (!locus || !Array.isArray(locus.genes)) throw new Error("every locus needs a genes array");
      for (const gene of locus.genes) {
        if (gene?.uid === undefined || genes.has(gene.uid)) {
          throw new Error("every gene needs a unique uid");
        }
        genes.add(gene.uid);
      }
    }
  }
  for (const group of data.groups) {
    if (!group || !Array.isArray(group.genes) || group.genes.some(uid => !genes.has(uid))) {
      throw new Error("every group must contain known gene IDs");
    }
  }
  for (const link of data.links) {
    if (!link?.query || !link?.target || !genes.has(link.query.uid) || !genes.has(link.target.uid)) {
      throw new Error("every link must refer to known genes");
    }
  }
}

export function parseProject(text) {
  let saved;
  try { saved = JSON.parse(text); }
  catch { throw new Error("the selected file is not valid JSON"); }
  if (saved?.format !== "clinker-project" || saved.version !== 2) throw new Error("the file is not a clinker project");
  if (!saved.config || typeof saved.config !== "object" || Array.isArray(saved.config)) {
    throw new Error("project appearance settings must be an object");
  }
  if (!saved.analysis) throw new Error("project is missing analysis data");
  validatePlotData(saved.data);
  return saved;
}

export function downloadText(filename, text, type) {
  const url = URL.createObjectURL(new Blob([text], { type }));
  const download = document.createElement("a");
  download.href = url;
  download.download = filename;
  document.body.append(download);
  download.click();
  download.remove();
  setTimeout(() => URL.revokeObjectURL(url), 0);
}
