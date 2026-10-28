import { mountPlotSelectionToolbar } from '../vendor/clustermap/editor.mjs?v=13c559062129a7eb';
import { createProject, validateProject } from './project.mjs?v=13c559062129a7eb';

const DEFAULT_CONFIG = {
  link: { bestOnly: true },
  plot: { renderer: 'webgpu' },
  legend: { columns: 4, position: 'bottom' },
};

// The host supplies elements and analysis results. This module never starts an
// analysis worker, reads a file, or decides which application controls to show.
export function createViewer({ plot, stage, editor, appearance = 'auto' }) {
  const makeChart = config => globalThis.ClusterMap.ClusterMap()
    .config(DEFAULT_CONFIG).config(config || {}).uiAppearance(appearance);
  let chart = makeChart();
  let analysis = null;
  let extensions = {};
  let disposeToolbar = null;
  let chartHost = null;
  editor.chart = chart;

  const editAnnotation = event => {
    if (event.detail?.chart !== chart || typeof editor.editAnnotations !== 'function') return;
    editor.open();
    editor.editAnnotations(event.detail.panelUid);
    event.preventDefault();
  };
  plot.addEventListener('clustermap:edit-annotation', editAnnotation);

  function configureMinimap(chart, keepEnabled) {
    const config = chart.config().plot;
    if (!['canvas', 'webgpu'].includes(config.renderer) || (keepEnabled && config.minimap.show)) return;
    const bounds = chart.figureBounds();
    if (!bounds) return;
    const width = bounds.maxX - bounds.minX;
    const height = bounds.maxY - bounds.minY;
    const viewport = plot.getBoundingClientRect();
    const camera = chart.state().camera;
    if (!viewport.width || !viewport.height) return;
    const left = bounds.minX * camera.k + camera.x;
    const top = bounds.minY * camera.k + camera.y;
    const show = left < -1 || top < -1
      || left + width * camera.k > viewport.width + 1
      || top + height * camera.k > viewport.height + 1;
    if (show !== config.minimap.show) {
      chart.config({ plot: { minimap: { show } } });
      chart.clearHistory();
    }
  }

  function mount(data, { config = chart.config(), state = null, analysis: session = null, keepMinimap = false } = {}) {
    const next = makeChart(config);
    const nextHost = plot.ownerDocument.createElement('div');
    nextHost.className = 'plot-chart';
    plot.append(nextHost);
    let nextToolbar = null;
    try {
      globalThis.d3.select(nextHost).datum(data).call(next);
      if (state) next.state(state);
      configureMinimap(next, keepMinimap);
      nextToolbar = mountPlotSelectionToolbar(next, stage);
      editor.chart = next;
    } catch (error) {
      nextToolbar?.();
      next.destroy();
      nextHost.remove();
      if (editor.chart !== chart) editor.chart = chart;
      throw error;
    }
    disposeToolbar?.();
    chart.destroy();
    chartHost?.remove();
    chart = next;
    chartHost = nextHost;
    disposeToolbar = nextToolbar;
    analysis = session;
    extensions = {};
    return chart;
  }

  return {
    get chart() { return chart; },
    get analysis() { return analysis; },
    get editorOpen() { return !editor.hidden; },
    mount,
    loadProject(saved) {
      validateProject(saved);
      mount(saved.data, { config: saved.config, state: saved.state, analysis: saved.analysis, keepMinimap: true });
      extensions = saved;
    },
    setAnalysis(session) { analysis = session; },
    setAppearance(preference) {
      appearance = preference;
      chart.uiAppearance(preference);
    },
    refreshEditor() { editor.chart = chart; },
    openEditor() { editor.open(); },
    closeEditor() { editor.close(); },
    openTreeAppearance() { return editor.openAppearance({ section: 'tree', focus: true }); },
    project(trees = {}) { return createProject(chart.project(), analysis, { ...extensions, ...trees }); },
    exportSvg() { return chart.exportSvg(); },
    destroy() {
      disposeToolbar?.();
      plot.removeEventListener('clustermap:edit-annotation', editAnnotation);
      editor.close();
      editor.chart = null;
      chart.destroy();
      chartHost?.remove();
    },
  };
}
