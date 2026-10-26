import { mountTreeInput } from './tree-input.mjs';
import { matchTreeLeaves, resolveTree, validateTree, identifyInputClusters } from './tree.mjs';

function frameTree(chart) {
  const bounds = chart.figureBounds();
  if (!bounds) return;
  const x = bounds.minX - 20;
  const y = bounds.minY - 20;
  const width = bounds.maxX - bounds.minX + 40;
  const height = bounds.maxY - bounds.minY + 40;
  const viewport = document.querySelector('#plot').getBoundingClientRect();
  const config = chart.config().plot;
  if (!(width > 0 && height > 0 && viewport.width > 0 && viewport.height > 0)) return;
  const fit = Math.min(viewport.width / width, viewport.height / height);
  const k = Math.min(config.maxZoom, Math.max(config.minZoom, config.renderer === 'svg' ? fit : Math.max(1, fit)));
  const camera = {
    x: Math.max(0, (viewport.width - width * k) / 2) - x * k,
    y: Math.max(0, (viewport.height - height * k) / 2) - y * k,
    k,
  };
  chart.state({ ...chart.state(), camera });
  chart.clearHistory();
}

export function mountTreeImport({ button, getChart, getRunInputs, onRunChange, onAppearance, onChange }) {
  const dialog = document.querySelector('#tree-dialog');
  const apply = dialog.querySelector('#apply-tree');
  const remove = document.querySelector('#remove-tree');
  const order = document.querySelector('#plot-order');
  const appearance = document.querySelector('#tree-appearance');
  const scaled = document.querySelector('#clustering-scaled');
  const scaleControl = document.querySelector('#clustering-scale-control');
  let inactive = null;
  let dendrogram = null;
  let activeSource = "imported";
  let boundChart;
  let unsubscribe;
  let target;
  let targetData;
  let inspectionVersion = 0;
  let continuation = null;
  const runInput = mountTreeInput(document.querySelector('#run-tree-fields'), {
    inline: true,
    onChange: () => onRunChange(runInput.ready()),
  });
  const plotInput = mountTreeInput(document.querySelector('#dialog-tree-fields'), {
    onChange: () => { apply.disabled = !plotInput.snapshot() || !plotInput.ready(); },
    isCurrent: () => !!continuation || (target === getChart() && targetData === target.data()),
  });
  const retainActiveTree = () => {
    const tree = getChart().data()?.tree;
    if (!tree) return;
    if (activeSource === 'clustering') dendrogram = structuredClone(tree);
    else inactive = structuredClone(tree);
  };
  const syncOrder = () => {
    const chart = getChart();
    const data = chart.data();
    if (!data) return;
    const clusterIds = new Set(data.clusters.map(cluster => cluster.uid));
    const compatible = tree => {
      if (!tree) return null;
      try { validateTree(tree, clusterIds); return tree; }
      catch { return null; }
    };
    inactive = compatible(inactive);
    dendrogram = compatible(activeSource === 'clustering' && data.tree ? data.tree : dendrogram);
    const displayed = chart.state().clusterOrder;
    const custom = displayed.length && displayed.some((uid, i) => uid !== data.clusters[i]?.uid);
    const selected = data.tree ? (activeSource === 'clustering' ? 'clustering' : 'tree') : custom ? 'manual' : 'computed';
    order.dataset.order = selected;
    scaleControl.hidden = selected !== 'clustering';
    scaled.checked = chart.config().tree.useBranchLengths !== false;
    scaled.disabled = !data.tree || !data.tree.nodes.filter(node => node.uid !== data.tree.root).every(node => Number.isFinite(node.branchLength));
    scaleControl.title = scaled.disabled ? 'Run a new analysis to obtain merge heights for this dendrogram.' : 'Use Ward linkage distances for horizontal spacing';
    const labels = { computed: 'Default', clustering: 'Similarity', tree: 'Tree', manual: 'Manual' };
    document.querySelector('#plot-order-label').textContent = labels[selected];
    for (const item of order.querySelectorAll('[data-order]')) item.setAttribute('aria-pressed', String(item.dataset.order === selected));
    const hasTree = (!!data.tree && activeSource !== 'clustering') || !!inactive;
    order.querySelector('[data-order="clustering"]').hidden = !dendrogram;
    order.querySelector('[data-order="tree"]').hidden = !hasTree;
    button.textContent = hasTree ? 'Replace tree (Newick)…' : 'Import tree (Newick)…';
    appearance.hidden = !data.tree && !inactive && !dendrogram;
    remove.hidden = !hasTree;
    remove.textContent = 'Remove imported tree';
    document.querySelector('#tree-controls-divider').hidden = appearance.hidden && remove.hidden;
  };
  const bind = () => {
    if (boundChart !== getChart()) {
      unsubscribe?.();
      inactive = null;
      dendrogram = null;
      activeSource = 'imported';
      boundChart = getChart();
      unsubscribe = boundChart.on('change', syncOrder);
    }
    syncOrder();
  };
  order.addEventListener('click', event => {
    const selected = event.target.closest('button[data-order]')?.dataset.order;
    if (!selected) return;
    order.open = false;
    const chart = getChart();
    try {
      if (selected === 'tree' || selected === 'clustering') {
        const next = selected === 'clustering' ? dendrogram : inactive;
        if (!next || (chart.data().tree && activeSource === (selected === 'clustering' ? 'clustering' : 'imported'))) return;
        validateTree(next, new Set(chart.data().clusters.map(cluster => cluster.uid)));
        retainActiveTree();
        activeSource = selected === 'clustering' ? 'clustering' : 'imported';
        chart.data({ ...chart.data(), tree: structuredClone(next) });
        if (selected === 'tree') inactive = null;
        frameTree(chart);
        onChange(selected === 'clustering' ? 'Showing clustering dendrogram from gene-cluster similarity.' : 'Using imported tree order.');
      } else if (selected === 'computed') {
        if (chart.data().tree) {
          retainActiveTree();
          const data = { ...chart.data() };
          delete data.tree;
          chart.data(data);
        }
        chart.state({ ...chart.state(), clusterOrder: [] });
        onChange('Using the computed default order. Trees are retained for switching back.');
      }
    } catch (error) { onChange(error.message || String(error)); }
    syncOrder();
  });
  scaled.addEventListener('change', () => {
    getChart().config({ tree: { useBranchLengths: scaled.checked } });
  });
  const inputsChanged = async () => {
    const version = ++inspectionVersion;
    runInput.setContext({ clusters: [], loading: true });
    try {
      const context = await getRunInputs();
      if (version === inspectionVersion) runInput.setContext(context);
    } catch (error) {
      if (version === inspectionVersion) runInput.setContext({ clusters: [], error: error.message });
    }
  };
  const open = (draft = null, clusters = null) => {
    target = getChart();
    targetData = target.data();
    if (!clusters && !targetData) return;
    plotInput.reset(draft);
    plotInput.setContext({ clusters: clusters || targetData.clusters });
    dialog.querySelector('#tree-title').textContent = !clusters && (targetData.tree || inactive) ? 'Replace tree' : 'Import tree';
    dialog.querySelector('#tree-import-note').hidden = !!clusters;
    dialog.showModal();
  };
  button.addEventListener('click', () => { order.open = false; open(); });
  apply.addEventListener('click', () => {
    if (!plotInput.ready()) return;
    const { tree, clusters, assignments } = plotInput.resolve();
    try {
      if (!continuation && (target !== getChart() || targetData !== target.data())) throw new Error('The figure changed. Reopen tree import before applying a tree.');
      const resolved = resolveTree(tree, clusters, assignments);
      if (continuation) {
        const waiting = continuation;
        continuation = null;
        runInput.reset(plotInput.snapshot());
        dialog.close();
        waiting.resolve(resolved);
      } else {
        retainActiveTree();
        activeSource = 'imported';
        target.data({ ...targetData, tree: resolved });
        inactive = null;
        syncOrder();
        frameTree(target);
        onChange('Attached tree. Click an internal node to rotate its children; row dragging follows the tree.');
        dialog.close();
      }
    } catch (error) { onChange(error.message || String(error)); }
  });
  remove.addEventListener('click', () => {
    order.open = false;
    const chart = getChart();
    if (activeSource !== 'clustering' && chart.data().tree) {
      const next = { ...chart.data() };
      delete next.tree;
      chart.data(next);
      chart.state({ ...chart.state(), clusterOrder: [] });
    }
    inactive = null;
    syncOrder();
    onChange('Removed imported tree.');
  });
  appearance.addEventListener('click', () => { order.open = false; onAppearance(); });
  dialog.addEventListener('keydown', event => event.stopPropagation());
  dialog.addEventListener('close', () => {
    if (continuation) { const waiting = continuation; continuation = null; waiting.workflow.cancel(); }
    plotInput.reset();
    const trigger = button.closest('details');
    if (!trigger?.closest('[inert]')) trigger?.querySelector('summary').focus();
  });
  runInput.reset();
  bind();
  return {
    close: () => { if (dialog.open) dialog.close(); },
    bind,
    projectTrees: () => ({
      inactiveTree: inactive ? structuredClone(inactive) : undefined,
      clusteringDendrogram: activeSource === 'clustering' && getChart().data()?.tree
        ? structuredClone(getChart().data().tree) : dendrogram ? structuredClone(dendrogram) : undefined,
      treeSource: getChart().data()?.tree && activeSource === 'clustering' ? 'clustering' : undefined,
    }),
    restoreTrees: saved => {
      inactive = saved.inactiveTree ? structuredClone(saved.inactiveTree) : null;
      dendrogram = saved.clusteringDendrogram ? structuredClone(saved.clusteringDendrogram) : null;
      activeSource = saved.treeSource || 'imported';
      syncOrder();
    },
    setDendrogram: tree => { dendrogram = tree || null; syncOrder(); },
    runTree: () => {
      if (!runInput.ready()) throw new Error('Resolve the tree input before starting analysis.');
      return runInput.snapshot();
    },
    runReady: runInput.ready,
    inputsChanged,
    resolveRunTree: async (clusters, draft, separateRecords, workflow) => {
      if (!draft) return null;
      const identified = identifyInputClusters(clusters, separateRecords);
      if (matchTreeLeaves(draft.tree, identified, draft.assignments).ready) return resolveTree(draft.tree, identified, draft.assignments);
      ++inspectionVersion;
      const pending = new Promise(resolve => {
        open(draft, identified);
        continuation = { resolve, workflow };
      });
      const abort = () => { continuation = null; if (dialog.open) dialog.close(); };
      workflow.signal.addEventListener('abort', abort, { once: true });
      try { return await workflow.wait(pending); }
      finally { workflow.signal.removeEventListener('abort', abort); }
    },
  };
}
