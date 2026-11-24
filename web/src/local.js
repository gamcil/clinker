import { createViewer } from './viewer/viewer.mjs';
import { mountViewerControls, readThemePreference } from './viewer/controls.mjs';
import { parseProject } from './viewer/project.mjs';
import { mountTreeImport } from './tree-import.mjs';
import { mountAnnotationImport } from './annotation-import.mjs';
import { createPlotStatus } from './viewer/status.mjs';
import { mountTour } from './viewer/tour.mjs';

const editor = document.querySelector('#plot-editor');
const viewer = createViewer({
  plot: document.querySelector('#plot'), stage: document.querySelector('.plot-stage'),
  editor, appearance: readThemePreference(),
});
const showStatus = createPlotStatus(document.querySelector('#plot-status'));
const tour = mountTour({ viewer, onError: showStatus });
let similarity = null;
let generation = 0;
const trees = mountTreeImport({
  button: document.querySelector('#manage-tree'), getChart: () => viewer.chart,
  onChange: (message, options) => { viewer.refreshEditor(); showStatus(message, options); },
});
mountAnnotationImport({
  button: document.querySelector('#import-annotations'), getChart: () => viewer.chart,
  onImported: showStatus,
});
mountViewerControls({
  viewer, editor, getTrees: () => trees.projectTrees(), getSimilarity: () => similarity,
  onError: showStatus,
});

async function nativeJson(path, options) {
  const response = await fetch(path, options);
  if (!response.ok) throw new Error(await response.text());
  return response.json();
}

function display(project, matrix) {
  trees.close();
  viewer.loadProject(project);
  trees.bind();
  trees.restoreTrees(project);
  viewer.closeEditor();
  similarity = matrix;
  document.querySelector('#plot-actions').hidden = false;
  document.querySelector('#edit-plot').hidden = false;
  for (const id of ['download-project', 'download-svg', 'download-matrix', 'edit-plot', 'manage-tree', 'import-annotations']) {
    document.getElementById(id).disabled = false;
  }
  showStatus(`Loaded figure · ${project.data.clusters.length.toLocaleString()} clusters · ${project.data.links.length.toLocaleString()} retained links`);
  tour.offer();
}

const fileInput = document.querySelector('#clinker-file');
document.querySelector('#load-clinker-file').addEventListener('click', () => fileInput.click());
fileInput.addEventListener('change', async () => {
  const file = fileInput.files[0];
  fileInput.value = '';
  if (!file) return;
  const current = ++generation;
  try {
    const project = parseProject(await file.text());
    const matrix = await nativeJson('api/similarity', {
      method: 'POST', headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(project.analysis),
    });
    if (current === generation) display(project, matrix);
  } catch (error) {
    if (current === generation) showStatus(`Could not open project: ${error.message}`, { tone: 'error' });
  }
});

try {
  const [project, matrix] = await Promise.all([nativeJson('api/project'), nativeJson('api/similarity')]);
  if (generation === 0) display(project, matrix);
} catch (error) {
  if (generation === 0) showStatus(`Could not open native results: ${error.message}`, { tone: 'error' });
}
