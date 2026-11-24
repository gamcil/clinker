import { driver } from '../vendor/driver/driver.js.mjs?v=f59ad68748ccfe73';
import { rangeShortcutKeys, formatShortcutKeys } from './shortcuts.mjs?v=f59ad68748ccfe73';
import { mountPlotSelectionToolbar } from '../vendor/clustermap/editor.mjs?v=f59ad68748ccfe73';

const seenKey = 'clinker-tour-seen';

export function mountTour({ viewer, onError }) {
  const prompt = document.querySelector('#tour-prompt');
  const help = document.querySelector('#help-dialog');
  const keys = rangeShortcutKeys(navigator.platform);
  const addKey = formatShortcutKeys(keys.add), removeKey = formatShortcutKeys(keys.remove);
  document.querySelector('#range-add-key').innerHTML = addKey;
  document.querySelector('#range-remove-key').innerHTML = removeKey;
  let seen = false;
  let starting = false;
  try { seen = localStorage.getItem(seenKey) === '1'; } catch {}
  function dismiss() {
    prompt.hidden = true;
    seen = true;
    try { localStorage.setItem(seenKey, '1'); } catch {}
  }
  async function start() {
    if (!viewer.chart.data() || starting) return;
    starting = true;
    let example;
    try {
      const response = await fetch(new URL('./tour-example.json', import.meta.url));
      if (!response.ok) throw new Error(`Could not load the tour example (${response.status})`);
      example = await response.json();
    } catch (error) {
      starting = false;
      onError(`Could not start the tour: ${error.message}`, { tone: 'error' });
      if (help.open) help.close();
      return;
    }
    dismiss();
    if (help.open) help.close();
    const editor = document.querySelector('#plot-editor');
    const originalChart = viewer.chart;
    const wasOpen = viewer.editorOpen;
    const previousTab = editor.shadowRoot?.querySelector('[data-editor-tabs]')?.value || 'data';
    editor.close();
    const stage = document.querySelector('.plot-stage');
    const plot = document.querySelector('#plot');
    const appBar = document.querySelector('.app-bar');
    const wasInert = appBar.inert;
    appBar.inert = true;
    const exampleHost = document.createElement('div');
    exampleHost.className = 'tour-example';
    const chartHost = document.createElement('div');
    chartHost.className = 'tour-example-chart';
    const label = document.createElement('div');
    label.className = 'tour-example-label';
    label.textContent = 'Example figure · changes here do not affect your project';
    exampleHost.append(chartHost, label);
    plot.append(exampleHost);
    stage.dataset.tourExample = '';
    let chart, disposeToolbar;
    try {
      chart = globalThis.ClusterMap.ClusterMap().config(originalChart.config())
        .config({ plot: { renderer: 'canvas', minimap: { show: false } }, tree: { show: true, scaleBar: { show: true } }, legend: { show: true } })
        .uiAppearance(document.documentElement.dataset.theme || 'auto');
      globalThis.d3.select(chartHost).datum(example).call(chart);
      editor.chart = chart;
      disposeToolbar = mountPlotSelectionToolbar(chart, stage, { editor });
      stage.lastElementChild.classList.add('tour-example-toolbar');
    } catch (error) {
      disposeToolbar?.(); chart?.destroy(); exampleHost.remove();
      delete stage.dataset.tourExample;
      appBar.inert = wasInert;
      editor.chart = originalChart;
      if (wasOpen) editor.open({ tab: previousTab });
      starting = false;
      onError(`Could not start the tour: ${error.message}`, { tone: 'error' });
      return;
    }
    starting = false;
    const exampleToolbar = stage.querySelector('.tour-example-toolbar');
    const exampleMenu = chartHost.querySelector('.tooltip');
    const plotControls = [exampleToolbar, exampleMenu];
    const positionPlotPopover = ({ wrapper }) => {
      if (!wrapper.classList.contains('tour-plot-popover')) return;
      const figure = chart.screenBounds({ figure: true });
      const target = document.querySelector('.tour-region-target, .tour-locus-target')?.getBoundingClientRect() || figure;
      const viewport = plot.getBoundingClientRect();
      if (!figure || !target || !viewport.width || !viewport.height) return;
      const { width, height } = wrapper.getBoundingClientRect();
      const margin = 12;
      const minX = Math.max(viewport.left, 0) + margin, minY = Math.max(viewport.top, 0) + margin;
      const maxX = Math.max(minX, Math.min(viewport.right, innerWidth) - width - margin);
      const maxY = Math.max(minY, Math.min(viewport.bottom, innerHeight) - height - margin);
      const clamp = (value, min, max) => Math.min(max, Math.max(min, value));
      const middleX = (target.left + target.right - width) / 2;
      const middleY = (target.top + target.bottom - height) / 2;
      const candidates = [
        [middleX, target.top - height - margin], [middleX, target.bottom + margin],
        [target.right + margin, middleY], [target.left - width - margin, middleY],
      ].map(([x, y]) => ({ x: clamp(x, minX, maxX), y: clamp(y, minY, maxY) }));
      const overlap = (p, rect) => Math.max(0, Math.min(p.x + width, rect.right) - Math.max(p.x, rect.left)) *
        Math.max(0, Math.min(p.y + height, rect.bottom) - Math.max(p.y, rect.top));
      const toolbar = !exampleToolbar.hidden ? exampleToolbar.getBoundingClientRect() : null;
      const score = p => overlap(p, target) * 4 + overlap(p, figure) + (toolbar ? overlap(p, toolbar) * 2 : 0);
      if (regionTarget?.dataset.region !== 'legend' || target.top - height - margin < minY) {
        candidates.sort((a, b) => score(a) - score(b));
      }
      wrapper.style.setProperty('left', `${candidates[0].x}px`, 'important');
      wrapper.style.setProperty('top', `${candidates[0].y}px`, 'important');
      wrapper.style.setProperty('transform', 'none', 'important');
    };
    let showingEditor = false;
    let locusTarget = null, regionTarget = null;
    const clearLocusTarget = () => {
      locusTarget?.remove();
      locusTarget = null;
    };
    const targetLocus = () => {
      if (locusTarget) return locusTarget;
      const loci = chart.visibleData()?.clusters.flatMap(cluster => cluster.loci.map(locus => locus.uid)) || [];
      if (!loci.length) return document.querySelector('#plot');
      const bounds = chart.screenBounds({ loci });
      if (!bounds || bounds.width <= 0 || bounds.height <= 0) {
        clearLocusTarget();
        return document.querySelector('#plot');
      }
      const viewport = document.querySelector('#plot').getBoundingClientRect();
      const left = Math.max(bounds.left, viewport.left), top = Math.max(bounds.top, viewport.top);
      const right = Math.min(bounds.right, viewport.right), bottom = Math.min(bounds.bottom, viewport.bottom);
      if (right <= left || bottom <= top) { clearLocusTarget(); return document.querySelector('#plot'); }
      locusTarget = document.createElement('div');
      locusTarget.className = 'tour-locus-target';
      Object.assign(locusTarget.style, { position: 'fixed', left: `${left}px`, top: `${top}px`, width: `${right - left}px`, height: `${bottom - top}px`, pointerEvents: 'none' });
      document.body.append(locusTarget);
      return locusTarget;
    };
    const clearRegionTarget = () => { regionTarget?.remove(); regionTarget = null; };
    const targetRegion = (targets, name, { matchLociHeight = false } = {}) => {
      if (regionTarget?.dataset.region === name) return regionTarget;
      clearRegionTarget();
      const bounds = chart.screenBounds(targets);
      if (!bounds) return plot;
      const viewport = plot.getBoundingClientRect();
      const loci = matchLociHeight ? chart.visibleData().clusters.flatMap(cluster => cluster.loci.map(locus => locus.uid)) : [];
      const vertical = loci.length ? chart.screenBounds({ loci }) : null;
      const left = Math.max(bounds.left, viewport.left), top = Math.max(Math.min(bounds.top, vertical?.top ?? bounds.top), viewport.top);
      const right = Math.min(bounds.right, viewport.right), bottom = Math.min(Math.max(bounds.bottom, vertical?.bottom ?? bounds.bottom), viewport.bottom);
      if (right <= left || bottom <= top) return plot;
      regionTarget = document.createElement('div');
      regionTarget.className = 'tour-region-target';
      regionTarget.dataset.region = name;
      Object.assign(regionTarget.style, { position: 'fixed', left: `${left}px`, top: `${top}px`, width: `${right - left}px`, height: `${bottom - top}px`, pointerEvents: 'none' });
      document.body.append(regionTarget);
      return regionTarget;
    };
    const showEditor = tab => {
      editor.open({ tab });
      showingEditor = true;
      return editor;
    };
    const restoreEditor = () => {
      if (!showingEditor) return;
      editor.close();
      showingEditor = false;
    };
    const steps = [
      ...(document.querySelector('#new-run')?.getClientRects().length ? [{ element: '#new-run', onHighlightStarted: () => { clearLocusTarget(); clearRegionTarget(); }, popover: { title: 'New run', description: 'Start a new comparison with your own files.' } }] : []),
      { element: '#file-menu', onHighlightStarted: () => { clearLocusTarget(); clearRegionTarget(); }, popover: { title: 'Save and export', description: 'Save or open a project, or export your figure and similarity table.' } },
      { element: '#add-data-menu', popover: { title: 'Add data', description: 'Import a tree or an annotation table to add to the figure.' } },
      { element: '#plot-order', onHighlightStarted: () => { clearRegionTarget(); restoreEditor(); }, popover: { title: 'Row order', description: 'Switch between the default order and available trees or similarity dendrogram.' } },
      { element: '#edit-plot', onHighlightStarted: restoreEditor, popover: { title: 'Figure editor', description: 'Use this button to open or close the figure editor.' } },
      { element: () => showEditor('data'), popover: { title: 'Edit figure data', description: 'Browse genes, groups, links and annotations. Select rows to edit several items together.' } },
      { element: () => showEditor('appearance'), popover: { title: 'Adjust appearance', description: 'Browse or search settings to change how the figure looks.' } },
      { element: '#appearance-menu', onHighlightStarted: restoreEditor, popover: { title: 'Theme', description: 'Choose the interface theme without changing the figure colours.' } },
      { element: '#plot-help', popover: { title: 'Help', description: 'Find shortcuts and reopen this tour here.' } },
      { element: '#plot', allowPageInteraction: true, allowedElements: plotControls, onHighlightStarted: () => { clearLocusTarget(); clearRegionTarget(); restoreEditor(); }, popover: { popoverClass: 'tour-plot-popover', title: 'Navigate', description: 'Scroll or pinch to zoom. Drag empty space to pan; hold <kbd>Space</kbd> to pan over a Canvas/GPU figure.' } },
      { element: targetLocus, allowPageInteraction: true, allowedElements: plotControls, popover: { popoverClass: 'tour-plot-popover', title: 'Select loci', description: `<kbd>Shift</kbd>-click a locus or gene to select or deselect it. From that locus, ${addKey} + <kbd>Shift</kbd>-click another to select the range; use ${removeKey} + <kbd>Shift</kbd>-click to deselect it.` } },
      { element: targetLocus, allowPageInteraction: true, allowedElements: plotControls, popover: { popoverClass: 'tour-plot-popover', title: 'Arrange loci', description: 'Drag a locus to move it, drag an end to trim it, or double-click it to flip its orientation. Use the selection toolbar for batch actions.' } },
      { element: () => {
        const ids = chart.visibleData().clusters.map(cluster => cluster.uid);
        return targetRegion({ clusterLabels: ids, clusterLocusLabels: ids }, 'labels', { matchLociHeight: true });
      }, allowPageInteraction: true, allowedElements: plotControls, onHighlightStarted: clearLocusTarget, popover: { popoverClass: 'tour-plot-popover', title: 'Arrange clusters', description: 'Drag a cluster label to reorder its row. Shift-click labels to select clusters and use the toolbar for batch actions.' } },
      { element: () => targetRegion({ tree: true, treeScaleBar: true }, 'tree', { matchLociHeight: true }), allowPageInteraction: true, allowedElements: plotControls, popover: { popoverClass: 'tour-plot-popover', title: 'Tree actions', description: 'Click a tree node to rotate branches, or right-click it for more actions. Choose a tree or similarity dendrogram under Row order.' } },
      { element: () => targetRegion({ legend: true }, 'legend'), allowPageInteraction: true, allowedElements: plotControls, popover: { popoverClass: 'tour-plot-popover', title: 'Legend', description: 'Click a group name or colour to edit it. Right-click for more group actions.' } },
    ];
    driver({ showProgress: true, allowClose: true, overlayOpacity: 0.15, steps, onPopoverPositioned: positionPlotPopover, onDestroyed: () => {
      clearLocusTarget(); clearRegionTarget(); disposeToolbar(); chart.destroy(); exampleHost.remove();
      delete stage.dataset.tourExample;
      appBar.inert = wasInert;
      editor.chart = originalChart;
      if (wasOpen) editor.open({ tab: previousTab });
      else editor.close();
    } }).drive();
  }
  prompt.querySelector('[data-tour-start]').addEventListener('click', start);
  prompt.querySelector('[data-tour-dismiss]').addEventListener('click', dismiss);
  document.querySelector('#help-tour').addEventListener('click', start);
  return { offer() {
    if (!viewer.chart.data()) return;
    document.querySelector('#help-tour').disabled = false;
    if (!seen) prompt.hidden = false;
  } };
}
