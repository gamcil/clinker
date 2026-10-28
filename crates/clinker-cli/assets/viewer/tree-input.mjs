import { Workflow } from './workflow.mjs?v=a2e7325c97ba540a';
import { parseNewick, matchTreeLeaves } from './viewer/tree.mjs?v=a2e7325c97ba540a';

// One set of fields serves the inline new-run input and the plot import dialog.
export function mountTreeInput(container, { inline = false, onChange = () => {}, isCurrent = () => true } = {}) {
  const fragment = document.querySelector('#tree-input-template').content.cloneNode(true);
  const prefix = inline ? 'run-' : '';
  for (const element of fragment.querySelectorAll('[id]')) element.id = prefix + element.id;
  container.append(fragment);
  if (inline) {
    container.querySelectorAll('.app-button, select.app-field').forEach(control => control.classList.add('compact-control'));
  }
  const find = id => container.querySelector(`#${prefix}${id}`);
  const input = find('tree-file');
  const choose = find('choose-tree-file');
  if (inline) choose.removeAttribute('autofocus');
  const fileName = find('tree-file-name');
  const clear = find('tree-clear');
  const status = find('tree-import-status');
  const matching = find('tree-matching');
  const rows = find('tree-matching-rows');
  matching.querySelector('th:last-child').textContent = inline ? 'Input cluster' : 'Plotted cluster';
  const missing = find('tree-missing');
  const rootMetadata = find('tree-root-metadata');
  const ignoreRoot = find('tree-ignore-root');
  const labels = find('tree-internal-labels');
  const scale = find('tree-support-scale');
  const options = find('tree-options');
  let context = { clusters: [] };
  let file = null;
  let text = '';
  let tree = null;
  let reading = null;
  let ready = true;
  let assignments = new Map();

  const report = (message = '', error = false) => {
    status.textContent = message;
    status.hidden = !message;
    status.dataset.tone = error ? 'error' : '';
  };
  const clusterLabel = cluster => {
    const name = cluster.source_file?.split(/[\\/]/).at(-1);
    const stem = name?.replace(/\.(gb|gbk|genbank|gbf|gbff|gff|gff3)$/i, '');
    return name && stem === cluster.name ? name : `${cluster.name || cluster.uid}${name ? ` · ${name}` : ''}`;
  };
  const renderMatches = () => {
    const result = matchTreeLeaves(tree, context.clusters, assignments);
    rows.replaceChildren();
    for (const row of result.rows) {
      const tr = document.createElement('tr');
      const leaf = document.createElement('td');
      leaf.textContent = row.label;
      const choice = document.createElement('td');
      const select = document.createElement('select');
      select.className = 'app-field flat-control compact-control';
      select.setAttribute('aria-label', `Cluster for tree leaf ${row.label}`);
      select.dataset.node = row.uid;
      select.append(new Option('Choose cluster…', ''));
      const chosen = context.clusters.find(cluster => cluster.uid === row.cluster);
      if (chosen) select.add(new Option(clusterLabel(chosen), chosen.uid));
      select.value = row.cluster;
      // Populate other choices on demand rather than storing N² options.
      const expand = () => {
        const selected = select.value;
        select.replaceChildren(new Option('Choose cluster…', ''));
        for (const cluster of context.clusters) select.add(new Option(clusterLabel(cluster), cluster.uid));
        select.value = selected;
        select.removeEventListener('focus', expand);
        select.removeEventListener('pointerdown', expand);
      };
      select.addEventListener('focus', expand);
      select.addEventListener('pointerdown', expand);
      select.addEventListener('change', () => {
        assignments.set(row.uid, select.value);
        refresh();
        rows.querySelector(`[data-node="${row.uid}"]`).focus();
      });
      choice.append(select);
      if (row.error) {
        const error = document.createElement('small');
        error.className = 'tree-match-error';
        error.textContent = row.error;
        choice.append(error);
      }
      tr.append(leaf, choice);
      rows.append(tr);
    }
    const matched = result.rows.filter(row => !row.error).length;
    matching.querySelector('summary').textContent = `${matched}/${result.rows.length} tree leaves matched`;
    missing.hidden = !result.missing.length;
    missing.textContent = result.missing.length ? `Clusters without a leaf: ${result.missing.map(clusterLabel).join('; ')}` : '';
    matching.hidden = false;
    ready = result.ready;
    if (!ready) {
      matching.open = true;
      if (!inline) report('Resolve the unmatched or duplicate leaves below.', true);
    }
  };
  const refresh = () => {
    ready = !file;
    tree = null;
    matching.hidden = true;
    rootMetadata.hidden = true;
    options.hidden = !file;
    clear.hidden = !inline || !file;
    choose.textContent = inline ? (file ? 'Change tree' : 'Choose tree') : 'Choose Newick file';
    fileName.textContent = file?.name || 'None selected';
    fileName.title = file?.name || '';
    scale.disabled = labels.value !== 'support';
    report();
    try {
      if (reading) { report(`Reading ${file.name}…`); return; }
      if (!file) return;
      if (!text.trim()) throw new Error('The tree file is empty.');
      const config = { internalLabels: labels.value, supportScale: scale.value };
      try { tree = parseNewick(text, config); }
      catch (error) {
        if (error.code !== 'ROOT_BRANCH_METADATA') throw error;
        rootMetadata.hidden = false;
        if (!ignoreRoot.checked) return;
        tree = parseNewick(text, { ...config, ignoreRootMetadata: true });
      }
      if (context.loading) report('Reading input files…');
      else if (context.error) report(context.error, true);
      else if (context.clusters.length) renderMatches();
      else if (inline) {
        ready = true;
        report(`${tree.nodes.filter(node => !node.children).length} tree leaves. Add input files to match them.`);
      }
    } catch (error) { tree = null; report(error.message || String(error), true); }
    finally { onChange(); }
  };
  const reset = (draft = null) => {
    reading?.dispose();
    reading = null;
    file = draft?.file || null;
    text = draft?.text || '';
    assignments = new Map(draft?.assignments);
    labels.value = draft?.options.internalLabels || 'support';
    scale.value = draft?.options.supportScale || 'unspecified';
    ignoreRoot.checked = draft?.options.ignoreRootMetadata || false;
    options.open = false;
    matching.open = !inline;
    input.value = '';
    refresh();
  };
  choose.addEventListener('click', () => input.click());
  clear.addEventListener('click', () => reset());
  input.addEventListener('change', async () => {
    const selected = input.files[0];
    reset();
    file = selected || null;
    if (!file) return;
    const workflow = new Workflow();
    reading = workflow;
    refresh();
    try {
      const contents = await workflow.readFile(file, true);
      workflow.throwIfCancelled();
      if (!isCurrent()) throw new Error('The figure changed. Reopen tree import to match the current clusters.');
      text = contents;
      reading = null;
      refresh();
    } catch (error) {
      if (!workflow.signal.aborted) { ready = false; report(error.message || String(error), true); onChange(); }
    } finally {
      workflow.dispose();
      if (reading === workflow) reading = null;
    }
  });
  for (const control of [labels, scale, ignoreRoot]) control.addEventListener('change', refresh);
  return {
    reset,
    setContext: value => { context = value; refresh(); },
    ready: () => ready && !reading,
    snapshot: () => file ? { file, text, tree: tree && structuredClone(tree), assignments: new Map(assignments), options: { internalLabels: labels.value, supportScale: scale.value, ignoreRootMetadata: ignoreRoot.checked } } : null,
    resolve: () => ({ tree, clusters: context.clusters, assignments }),
  };
}
