import { mountTreeInput } from './tree-input.mjs?v=20261018-tree-v8';
import { matchTreeLeaves, resolveTree, validateTree, identifyInputClusters } from './tree.mjs?v=20261018-tree-v5';

function frameTree(chart) {
  // Use public export bounds so framing works across SVG, Canvas and WebGPU.
  const svg = new DOMParser().parseFromString(chart.exportSvg(), 'image/svg+xml').documentElement;
  const [x, y, width, height] = svg.getAttribute('viewBox').split(/\s+/).map(Number);
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
  let inactive = null;
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
  const syncOrder = () => {
    const chart = getChart();
    const data = chart.data();
    if (!data) return;
    const displayed = chart.state().clusterOrder;
    const custom = displayed.length && displayed.some((uid, i) => uid !== data.clusters[i]?.uid);
    const selected = data.tree ? 'tree' : custom ? 'manual' : 'computed';
    order.dataset.order = selected;
    document.querySelector('#plot-order-label').textContent = `Order: ${selected[0].toUpperCase() + selected.slice(1)}`;
    for (const item of order.querySelectorAll('[data-order]')) item.setAttribute('aria-pressed', String(item.dataset.order === selected));
    const hasTree = !!data.tree || !!inactive;
    order.querySelector('[data-order="tree"]').hidden = !hasTree;
    button.textContent = hasTree ? 'Replace tree…' : 'Import tree…';
    appearance.hidden = !hasTree;
    remove.hidden = !hasTree;
  };
  const bind = () => {
    if (boundChart !== getChart()) {
      unsubscribe?.();
      inactive = null;
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
      if (selected === 'tree') {
        if (!inactive) return;
        validateTree(inactive, new Set(chart.data().clusters.map(cluster => cluster.uid)));
        const tree = inactive;
        chart.data({ ...chart.data(), tree });
        inactive = null;
        frameTree(chart);
        onChange('Using tree order.');
      } else if (selected === 'computed') {
        if (chart.data().tree) {
          inactive = structuredClone(chart.data().tree);
          const data = { ...chart.data() };
          delete data.tree;
          chart.data(data);
        }
        chart.state({ ...chart.state(), clusterOrder: [] });
        onChange('Using the computed default order. Any attached tree is retained for switching back.');
      }
    } catch (error) { onChange(error.message || String(error)); }
    syncOrder();
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
    const next = { ...chart.data() };
    delete next.tree;
    chart.data(next);
    chart.state({ ...chart.state(), clusterOrder: [] });
    inactive = null;
    syncOrder();
    onChange('Removed tree. Using the computed default order.');
  });
  appearance.addEventListener('click', () => { order.open = false; onAppearance(); });
  dialog.addEventListener('keydown', event => event.stopPropagation());
  dialog.addEventListener('close', () => {
    if (continuation) { const waiting = continuation; continuation = null; waiting.workflow.cancel(); }
    plotInput.reset();
    order.querySelector('summary').focus();
  });
  runInput.reset();
  bind();
  return {
    close: () => { if (dialog.open) dialog.close(); },
    bind,
    inactiveTree: () => inactive ? structuredClone(inactive) : undefined,
    restoreInactiveTree: tree => { inactive = tree ? structuredClone(tree) : null; syncOrder(); },
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
