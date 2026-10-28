import { mountPlotSelectionToolbar } from '../vendor/clustermap/editor.mjs?v=65753b40bd6366c3';
import { createProject, validateProject } from './project.mjs?v=65753b40bd6366c3';

const DEFAULT_CONFIG = {
  link: { bestOnly: true },
  plot: { renderer: 'webgpu' },
  legend: { columns: 4, position: 'bottom' },
};

// The host supplies elements and analysis results. This module never starts an
// analysis worker, reads a file, or decides which application controls to show.
export function createViewer({ plot, stage, editor, appearance = 'auto' }) {
  const selection = globalThis.d3.select(plot);
  const makeChart = config => globalThis.ClusterMap.ClusterMap()
    .config(DEFAULT_CONFIG).config(config || {}).uiAppearance(appearance);
  let chart = makeChart();
  let analysis = null;
  let extensions = {};
  let disposeToolbar = null;
  editor.chart = chart;

  const editAnnotation = event => {
    if (event.detail?.chart !== chart || typeof editor.editAnnotations !== 'function') return;
    editor.open();
    editor.editAnnotations(event.detail.panelUid);
    event.preventDefault();
  };
  plot.addEventListener('clustermap:edit-annotation', editAnnotation);

  function configureMinimap(keepEnabled) {
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
    disposeToolbar?.();
    disposeToolbar = null;
    const previous = chart;
    chart = next;
    previous.destroy();
    analysis = session;
    extensions = {};
    selection.datum(data).call(chart);
    if (state) chart.state(state);
    configureMinimap(keepMinimap);
    editor.chart = chart;
    disposeToolbar = mountPlotSelectionToolbar(chart, stage);
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
    },
  };
}
