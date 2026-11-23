import { driver } from '../vendor/driver/driver.js.mjs?v=90af91a815e4c2ca';
import { rangeShortcutKeys, formatShortcutKeys } from './shortcuts.mjs?v=90af91a815e4c2ca';

const seenKey = 'clinker-tour-seen';

export function mountTour({ viewer }) {
  const prompt = document.querySelector('#tour-prompt');
  const help = document.querySelector('#help-dialog');
  const keys = rangeShortcutKeys(navigator.platform);
  const addKey = formatShortcutKeys(keys.add), removeKey = formatShortcutKeys(keys.remove);
  document.querySelector('#range-add-key').innerHTML = addKey;
  document.querySelector('#range-remove-key').innerHTML = removeKey;
  let seen = false;
  try { seen = localStorage.getItem(seenKey) === '1'; } catch {}
  function dismiss() {
    prompt.hidden = true;
    seen = true;
    try { localStorage.setItem(seenKey, '1'); } catch {}
  }
  function start() {
    if (!viewer.chart.data()) return;
    dismiss();
    if (help.open) help.close();
    const editor = document.querySelector('#plot-editor');
    const wasOpen = viewer.editorOpen;
    const previousTab = editor.shadowRoot?.querySelector('[data-editor-tabs]')?.value || 'data';
    let showingEditor = false;
    const showEditor = tab => {
      editor.open({ tab });
      showingEditor = true;
      return editor;
    };
    const restoreEditor = () => {
      if (!showingEditor) return;
      editor.open({ tab: previousTab });
      if (!wasOpen) editor.close();
      showingEditor = false;
    };
    const steps = [
      { element: '#plot', popover: { popoverClass: 'tour-plot-popover', title: 'Navigate', description: 'Scroll or pinch to zoom. Drag empty space to pan.' } },
      { element: '#plot', popover: { popoverClass: 'tour-plot-popover', title: 'Select loci', description: `<kbd>Shift</kbd>-click a locus or gene to select or deselect it. From that locus, ${addKey} + <kbd>Shift</kbd>-click another to select the range; use ${removeKey} + <kbd>Shift</kbd>-click to deselect it.` } },
      ...(document.querySelector('#new-run')?.getClientRects().length ? [{ element: '#new-run', popover: { title: 'New run', description: 'Start a new comparison with your own files.' } }] : []),
      { element: '#file-menu', popover: { title: 'Save and export', description: 'Save or open a project, or export your figure and similarity table.' } },
      { element: '#add-data-menu', popover: { title: 'Add data', description: 'Import a tree or an annotation table to add to the figure.' } },
      { element: '#plot-order', onHighlightStarted: restoreEditor, popover: { title: 'Row order', description: 'Switch between the default order and available trees or similarity dendrogram.' } },
      { element: '#edit-plot', onHighlightStarted: restoreEditor, popover: { title: 'Figure editor', description: 'Use this button to open or close the figure editor.' } },
      { element: () => showEditor('data'), popover: { title: 'Edit figure data', description: 'Browse genes, groups, links and annotations. Select rows to edit several items together.' } },
      { element: () => showEditor('appearance'), popover: { title: 'Adjust appearance', description: 'Browse or search settings to change how the figure looks.' } },
      { element: '#appearance-menu', onHighlightStarted: restoreEditor, popover: { title: 'Theme', description: 'Choose the interface theme without changing the figure colours.' } },
      { element: '#plot-help', popover: { title: 'Help', description: 'Find shortcuts and reopen this tour here.' } },
    ];
    driver({ showProgress: true, allowClose: true, steps, onDestroyed: restoreEditor }).drive();
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
