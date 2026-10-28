import { validateTree } from './tree.mjs?v=cc18ab0a8d958d5c';

export const PROJECT_FORMAT = 'clinker-project';
export const PROJECT_VERSION = 2;

export function displayedClusterOrder({ data, state }) {
  if (data.tree !== undefined) return validateTree(data.tree, new Set(data.clusters.map(cluster => cluster.uid)));
  return state?.clusterOrder?.length ? state.clusterOrder : data.clusters.map(cluster => cluster.uid);
}

function csvCell(value) {
  const text = String(value);
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text;
}

export function similarityCsv(similarityMatrix, project) {
  const clusters = new Map(project.data.clusters.map(cluster => [cluster.uid, cluster]));
  const order = displayedClusterOrder(project);
  const displayed = order.map(uid => {
    const index = /^cluster-(0|[1-9]\d*)$/.test(uid) ? Number(uid.slice(8)) : -1;
    if (!clusters.has(uid) || index < 0 || index >= similarityMatrix.length) {
      throw new Error("Figure and analysis cluster IDs do not match.");
    }
    return { index, name: clusters.get(uid).name };
  });
  const plotPosition = new Map(displayed.map(({ index }, position) => [index, position + 1]));
  displayed.sort((a, b) => a.index - b.index);
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
  for (let i = 0; i < displayed.length; i += 1) {
    const query = displayed[i];
    for (let j = i + 1; j < displayed.length; j += 1) {
      const target = displayed[j];
      const pair = similarityMatrix[query.index][target.index];
      rows.push([
        plotPosition.get(query.index),
        query.name,
        query.index,
        plotPosition.get(target.index),
        target.name,
        target.index,
        pair.queryCoverage,
        pair.targetCoverage,
        pair.similarity,
        1 - pair.similarity,
      ]);
    }
  }
  return `${rows.map(row => row.map(csvCell).join(",")).join("\n")}\n`;
}

function uniqueId(record, ids, kind) {
  if (typeof record?.uid !== "string" || !record.uid || ids.has(record.uid)) {
    throw new Error(`every ${kind} needs a unique uid`);
  }
  ids.add(record.uid);
}

function validateSpan(record) {
  if (!Number.isFinite(record.start) || !Number.isFinite(record.end) || record.end < record.start) {
    throw new Error("sequence and gene coordinates must be finite and ordered");
  }
}

function validatePlotData(data, analysis) {
  if (!data || typeof data !== "object" || Array.isArray(data)) {
    throw new Error("the file must contain a plot-data object");
  }
  if (!Array.isArray(data.clusters) || !Array.isArray(data.links) || !Array.isArray(data.groups)) {
    throw new Error("plot data must contain clusters, links, and groups arrays");
  }
  const clusters = new Set();
  const loci = new Set();
  const genes = new Set();
  const groups = new Set();
  const links = new Set();
  const sourceClusters = analysis?.analysis?.clusters;
  if (!Array.isArray(sourceClusters)) throw new Error("project is missing analysis clusters");
  const knownGenes = new Set();
  sourceClusters.forEach((cluster, c) => {
    if (!Array.isArray(cluster?.loci)) throw new Error("analysis clusters need loci arrays");
    cluster.loci.forEach((locus, l) => {
      if (!Array.isArray(locus?.genes)) throw new Error("analysis loci need genes arrays");
      locus.genes.forEach((_, g) => knownGenes.add(`gene-${c}-${l}-${g}`));
    });
  });
  for (const cluster of data.clusters) {
    if (!cluster || !Array.isArray(cluster.loci)) throw new Error("every cluster needs a loci array");
    uniqueId(cluster, clusters, "cluster");
    const index = /^cluster-(0|[1-9]\d*)$/.test(cluster.uid) ? Number(cluster.uid.slice(8)) : -1;
    const source = sourceClusters[index];
    if (!source || !Array.isArray(source.loci)) throw new Error("figure and analysis clusters do not match");
    const sourceLoci = new Map(source.loci.map((locus, i) => [`locus-${index}-${i}`, locus]));
    for (const locus of cluster.loci) {
      if (!locus || !Array.isArray(locus.genes)) throw new Error("every locus needs a genes array");
      uniqueId(locus, loci, "locus");
      validateSpan(locus);
      if (locus.bio) validateSpan(locus.bio);
      const sourceLocus = sourceLoci.get(locus.uid);
      if (!Array.isArray(sourceLocus?.genes)) throw new Error("figure and analysis loci do not match");
      const sourceGenes = new Set(sourceLocus.genes.map((_, i) => `gene-${locus.uid.slice(6)}-${i}`));
      for (const gene of locus.genes) {
        uniqueId(gene, genes, "gene");
        validateSpan(gene);
        if (gene.bio) validateSpan(gene.bio);
        if (!sourceGenes.has(gene.uid)) throw new Error("figure and analysis genes do not match");
      }
    }
  }
  for (const group of data.groups) {
    uniqueId(group, groups, "group");
    if (!group || !Array.isArray(group.genes) || group.genes.some(uid => !knownGenes.has(uid))) {
      throw new Error("every group must contain known gene IDs");
    }
  }
  for (const link of data.links) {
    uniqueId(link, links, "link");
    // The chart retains links to deleted genes for editing/undo. They still
    // need valid analysis IDs, but need not have visible endpoints.
    if (!link?.query || !link?.target || !knownGenes.has(link.query.uid) || !knownGenes.has(link.target.uid)) {
      throw new Error("every link must refer to known genes");
    }
  }
  return clusters;
}

export function parseProject(text) {
  let saved;
  try { saved = JSON.parse(text); }
  catch { throw new Error("the selected file is not valid JSON"); }
  return validateProject(saved);
}

export function validateProject(saved) {
  if (saved?.format !== PROJECT_FORMAT || saved.version !== PROJECT_VERSION) throw new Error("the file is not a clinker project");
  if (!saved.config || typeof saved.config !== "object" || Array.isArray(saved.config)) {
    throw new Error("project appearance settings must be an object");
  }
  if (!saved.analysis) throw new Error("project is missing analysis data");
  const clusterIds = validatePlotData(saved.data, saved.analysis);
  const treeOrder = saved.data.tree === undefined ? null : validateTree(saved.data.tree, clusterIds);
  if (saved.inactiveTree !== undefined) {
    if (treeOrder && saved.treeSource !== "clustering") throw new Error("project cannot have both an active and inactive tree");
    validateTree(saved.inactiveTree, clusterIds);
  }
  if (saved.clusteringDendrogram !== undefined) validateTree(saved.clusteringDendrogram, clusterIds);
  if (saved.treeSource !== undefined && (saved.treeSource !== 'clustering' || !treeOrder || !saved.clusteringDendrogram)) {
    throw new Error('Invalid clustering tree source');
  }
  if (saved.state != null) {
    if (typeof saved.state !== "object" || Array.isArray(saved.state) || saved.state.version !== 1) {
      throw new Error("unsupported project layout state");
    }
    const order = saved.state.clusterOrder;
    if (!Array.isArray(order) || (order.length > 0
        && (order.length !== clusterIds.size || new Set(order).size !== order.length
          || order.some(uid => !clusterIds.has(uid))))) {
      throw new Error("project cluster order must contain each displayed cluster once");
    }
    if (treeOrder && order.length && order.some((uid, i) => uid !== treeOrder[i])) {
      throw new Error("project cluster order must agree with its tree child order");
    }
  }
  return saved;
}

// Preserve extension fields from an opened project; the current figure and
// retained analysis replace their saved counterparts. UI appearance is separate.
export function createProject(figure, analysis, extensions = {}) {
  return validateProject({ ...extensions, ...figure, format: PROJECT_FORMAT, version: PROJECT_VERSION, analysis });
}

export function serializeProject(project) {
  return `${JSON.stringify(validateProject(project), null, 2)}\n`;
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
