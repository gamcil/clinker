import test from 'node:test';
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const root = resolve(process.env.CLINKER_WEB_ROOT || fileURLToPath(new URL('../../', import.meta.url)));
const browserName = process.env.CLINKER_TEST_BROWSER || 'chromium';
const browserType = { chromium, firefox, webkit }[browserName];
if (!browserType) throw new Error(`Unknown test browser: ${browserName}`);

test('analysis, cancellation, failure recovery and project exports preserve the displayed figure', { timeout: 120000 }, async t => {
  const server = createServer(async (request, response) => {
    const path = resolve(root, '.' + decodeURIComponent(new URL(request.url, 'http://localhost').pathname).replace(/\/$/, '/index.html'));
    if (!path.startsWith(root + '/')) { response.writeHead(403).end(); return; }
    try {
      const bytes = await readFile(path);
      response.setHeader('Content-Type', { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm' }[extname(path)] || 'application/octet-stream');
      response.end(bytes);
    } catch { response.writeHead(404).end(); }
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  t.after(() => { server.closeAllConnections(); return new Promise(resolve => server.close(resolve)); });
  const browser = await browserType.launch({ headless: true });
  t.after(() => browser.close());
  const page = await browser.newPage();
  await page.addInitScript(() => {
    window.inputAudit = { reads: [], aborts: 0 };
    const readAsText = FileReader.prototype.readAsText;
    FileReader.prototype.readAsText = function(file, ...args) {
      window.inputAudit.reads.push(file.name);
      this.addEventListener('abort', () => { window.inputAudit.aborts += 1; }, { once: true });
      return readAsText.call(this, file, ...args);
    };
    const NativeWorker = window.Worker;
    window.workerAudit = { hold: null, pending: null, live: new Set() };
    window.Worker = class extends NativeWorker {
      constructor(...args) {
        super(...args);
        window.workerAudit.live.add(this);
      }
      postMessage(message, ...args) {
        if (message.type === window.workerAudit.hold) {
          window.workerAudit.pending = { worker: this, lateMessage: this.onmessage };
          return;
        }
        super.postMessage(message, ...args);
      }
      terminate() {
        window.workerAudit.live.delete(this);
        super.terminate();
      }
    };
  });
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  const helpButton = page.getByRole('button', { name: 'Plot help', exact: true });
  const help = page.getByRole('dialog', { name: 'Plot help' });
  await helpButton.click();
  assert.equal(await help.isVisible(), true);
  assert.equal(await help.getByRole('heading', { name: 'Select loci' }).isVisible(), true);
  assert.equal(await help.getByRole('button', { name: 'Close help' }).evaluate(el => el === document.activeElement), true);
  await page.keyboard.press('Tab');
  assert.equal(await help.getByRole('button', { name: 'Move around', exact: true }).evaluate(el => el === document.activeElement), true);
  await page.keyboard.press('Shift+Tab');
  assert.equal(await help.getByRole('button', { name: 'Close help' }).evaluate(el => el === document.activeElement), true);
  await page.keyboard.press('Escape');
  assert.equal(await help.isVisible(), false);
  assert.equal(await helpButton.evaluate(el => el === document.activeElement), true);
  await page.setViewportSize({ width: 390, height: 844 });
  await helpButton.click();
  const dimensions = await help.getByRole('region', { name: 'Plot controls and shortcuts' }).evaluate(el => ({
    scrolls: el.scrollHeight > el.clientHeight,
    fits: el.scrollWidth <= el.clientWidth,
  }));
  assert.deepEqual(dimensions, { scrolls: true, fits: true });
  await help.getByRole('button', { name: 'Save your work', exact: true }).click();
  assert.equal(await help.getByRole('region', { name: 'Plot controls and shortcuts' }).evaluate(el => el.scrollTop > 0), true);
  await page.mouse.click(2, 2);
  assert.equal(await help.isVisible(), false);
  await page.setViewportSize({ width: 1280, height: 720 });
  await page.evaluate(() => {
    const transfer = new DataTransfer();
    transfer.items.add(new File(['LOCUS       old 10 bp DNA\n', ' '.repeat(1024 * 1024)], 'old.gbk'));
    transfer.items.add(new File(['LOCUS       queued 10 bp DNA'], 'queued.gbk'));
    const input = document.querySelector('#files');
    input.files = transfer.files;
    input.dispatchEvent(new Event('change'));
    document.querySelector('#clear-input-files').click();
  });
  await page.waitForFunction(() => !document.querySelector('#record-ranges').hasAttribute('aria-busy'));
  assert.deepEqual(await page.evaluate(() => window.inputAudit), { reads: ['old.gbk'], aborts: 1 });
  assert.equal(await page.locator('#analyse').isDisabled(), true);
  assert.equal(await page.locator('[data-record-range]').count(), 0);
  await page.locator('#files').setInputFiles([
    resolve(root, 'examples/A. versicolor CBS 583.65.gbk'),
    resolve(root, 'examples/A. burnettii MST-FP2249.gbk'),
    resolve(root, 'examples/P. vexata CBS 129021.gbk'),
  ]);
  await page.waitForFunction(() => !document.querySelector('#analyse').disabled);
  const inspectionReads = await page.evaluate(() => window.inputAudit.reads);
  await page.locator('#files').dispatchEvent('change');
  await page.waitForFunction(() => !document.querySelector('#analyse').disabled);
  assert.deepEqual(await page.evaluate(() => window.inputAudit.reads), inspectionReads);
  await page.locator('#analyse').click();
  await page.waitForFunction(() => document.querySelector('#run-panel').hidden, {}, { timeout: 60000 });
  const expected = await page.evaluate(() => {
    const chart = document.querySelector('#plot-editor').chart;
    const gene = chart.data().clusters[0].loci[0].genes[0];
    chart.patch([
      { type: 'genes.update', ids: [gene.uid], changes: { label: 'Audit edited gene' } },
      { type: 'clusters.update', ids: [chart.data().clusters[0].uid], changes: { name: 'Audit renamed cluster' } },
    ]);
    chart.config({ plot: { renderer: 'svg' }, gene: { label: { show: true } } });
    const state = chart.state();
    state.clusterOrder = chart.data().clusters.map(c => c.uid).reverse();
    chart.state(state);
    return { gene: gene.uid, state: chart.state(), names: Object.fromEntries(chart.data().clusters.map(c => [c.uid, c.name])) };
  });
  const figure = () => page.evaluate(() => document.querySelector('#plot-editor').chart.project());
  const beforeDelete = await figure();
  const selectedLocus = beforeDelete.data.clusters[0].loci[0].uid;
  await page.evaluate(uid => document.querySelector('#plot-editor').chart.locusSelection([uid]), selectedLocus);
  await helpButton.click();
  await page.keyboard.press('Escape');
  assert.deepEqual(await page.evaluate(() => document.querySelector('#plot-editor').chart.locusSelection()), [selectedLocus]);
  await page.evaluate(() => document.querySelector('#plot-editor').chart.locusSelection([]));
  await page.locator('#edit-plot').click();
  const editor = page.locator('#plot-editor');
  await editor.locator('clinker-view-switcher[data-kind-picker]').getByRole('button', { name: 'Genes', exact: true }).click();
  const targetCluster = beforeDelete.data.clusters.find(cluster => cluster.loci.length === 1 && cluster.loci[0].genes.length);
  assert.ok(targetCluster, 'the example has a single-locus cluster for bulk gene selection');
  const deletedGenes = targetCluster.loci[0].genes.map(gene => gene.uid);
  const selection = editor.getByRole('checkbox', {
    name: `Select ${deletedGenes.length} genes in ${targetCluster.label || targetCluster.name || targetCluster.uid}`, exact: true,
  });
  assert.ok(deletedGenes.length > 0);
  await selection.check();
  await editor.locator('clinker-record-actions').getByText('Remove…', { exact: true }).click();
  page.once('dialog', dialog => dialog.accept());
  await editor.getByRole('button', { name: 'Delete selected genes', exact: true }).click();
  const afterDelete = await figure();
  const geneIds = snapshot => snapshot.data.clusters.flatMap(c => c.loci.flatMap(l => l.genes.map(g => g.uid)));
  const contents = ({ data }) => ({
    clusters: data.clusters, links: data.links, groups: data.groups,
    automaticGroups: data.config?.updateGroups !== false,
  });
  assert.deepEqual(geneIds(afterDelete), geneIds(beforeDelete).filter(uid => !deletedGenes.includes(uid)));
  await page.locator('#plot-editor [data-command="undo"]').click();
  assert.deepEqual(contents(await figure()), contents(beforeDelete));
  await page.locator('#plot-editor [data-command="redo"]').click();
  const afterRedo = await figure();
  assert.deepEqual(contents(afterRedo), contents(afterDelete));
  assert.deepEqual(afterRedo.state.clusterOrder, expected.state.clusterOrder);
  await page.locator('#edit-plot').click();
  async function download(id) {
    await page.locator('#file-menu > summary').click();
    const downloading = page.waitForEvent('download');
    await page.locator(id).click();
    return readFile(await (await downloading).path(), 'utf8');
  }
  function checkCsv(csv) {
    for (const line of csv.trim().split('\n').slice(1)) {
      const row = line.split(',');
      assert.equal(Number(row[0]), expected.state.clusterOrder.indexOf(`cluster-${row[2]}`) + 1);
      assert.equal(Number(row[3]), expected.state.clusterOrder.indexOf(`cluster-${row[5]}`) + 1);
      assert.equal(row[1], expected.names[`cluster-${row[2]}`]);
      assert.equal(row[4], expected.names[`cluster-${row[5]}`]);
    }
  }
  checkCsv(await download('#download-matrix'));
  const project = await download('#download-project');
  const openProject = buffer => page.locator('#clinker-file').setInputFiles({
    name: 'edited.json', mimeType: 'application/json', buffer: Buffer.from(buffer),
  });
  async function checkRecovery() {
    await page.waitForFunction(() => !document.querySelector('.run-card').classList.contains('is-busy'));
    assert.equal(await page.evaluate(() => window.workerAudit.live.size), 0);
    for (const id of ['#analyse', '#load-example', '#open-project-start', '#download-matrix', '#download-svg', '#download-project']) {
      assert.equal(await page.locator(id).isEnabled(), true, id);
    }
    const preserved = await page.evaluate(() => document.querySelector('#plot-editor').chart.project());
    assert.deepEqual(preserved.data, JSON.parse(project).data);
    assert.deepEqual(preserved.state.clusterOrder, expected.state.clusterOrder);
  }
  await page.locator('#new-run').click();
  // Hold each worker stage at its request boundary so cancellation is reliable
  // even when the example analysis finishes in milliseconds.
  for (const stage of ['parse-start', 'parse-finish', 'analyse-batch', 'post-process', 'open-session']) {
    await page.evaluate(stage => {
      window.workerAudit.hold = stage;
      window.workerAudit.pending = null;
    }, stage);
    if (stage === 'open-session') await openProject(project);
    else await page.locator('#analyse').click();
    await page.waitForFunction(() => window.workerAudit.pending !== null);
    await page.locator('#cancel-run').click();
    await checkRecovery();
    await page.evaluate(() => window.workerAudit.pending.lateMessage({
      data: { type: 'error', message: 'Abandoned worker response' },
    }));
    assert.match(await page.locator('#run-status').textContent(), /Cancelled/);
  }
  for (const stage of ['parse-start', 'analyse-batch', 'post-process', 'open-session']) {
    await page.evaluate(stage => {
      window.workerAudit.hold = stage;
      window.workerAudit.pending = null;
    }, stage);
    if (stage === 'open-session') await openProject(project);
    else await page.locator('#analyse').click();
    await page.waitForFunction(() => window.workerAudit.pending !== null);
    await page.evaluate(() => window.workerAudit.pending.worker.onmessageerror());
    await checkRecovery();
    assert.match(await page.locator('#run-status').textContent(), /Could not read/);
  }
  await page.evaluate(() => { window.workerAudit.hold = null; });
  await openProject('{invalid json');
  await checkRecovery();
  assert.match(await page.locator('#run-status').textContent(), /Could not open project/);
  const invalid = JSON.parse(project);
  invalid.data.clusters[1].uid = invalid.data.clusters[0].uid;
  await openProject(JSON.stringify(invalid));
  await checkRecovery();
  assert.match(await page.locator('#run-status').textContent(), /unique uid/);
  const incompatible = JSON.parse(project);
  incompatible.analysis.version = 999;
  await openProject(JSON.stringify(incompatible));
  await checkRecovery();
  assert.match(await page.locator('#run-status').textContent(), /unsupported analysis-session version/);
  await openProject(project);
  await page.waitForFunction(() => document.querySelector('#run-panel').hidden, {}, { timeout: 60000 });
  await page.reload();
  await page.locator('#clinker-file').setInputFiles({ name: 'edited.json', mimeType: 'application/json', buffer: Buffer.from(project) });
  await page.waitForFunction(() => document.querySelector('#plot-status').textContent.includes('Loaded clinker project'), {}, { timeout: 60000 });
  const reopened = await page.evaluate(() => document.querySelector('#plot-editor').chart.project());
  assert.deepEqual(reopened.data, afterRedo.data);
  assert.deepEqual(reopened.state.clusterOrder, expected.state.clusterOrder);
  assert.equal(reopened.data.clusters.flatMap(c => c.loci.flatMap(l => l.genes)).find(g => g.uid === expected.gene).label, 'Audit edited gene');
  checkCsv(await download('#download-matrix'));
  assert.match(await download('#download-svg'), /Audit edited gene/);
  await page.setViewportSize({ width: 2000, height: 1400 });
  await page.reload();
  await page.locator('#load-example').click();
  await page.waitForFunction(() => !document.querySelector('#edit-plot').disabled);
  assert.equal(await page.evaluate(() => document.querySelector('#plot-editor').chart.config().plot.minimap.show), false, 'a fitting plot needs no minimap');
  await page.locator('#edit-plot').click();
  assert.ok((await page.locator('#plot-editor').boundingBox()).height >= 490, 'the editor opens with a useful working height');
  await page.locator('#edit-plot').click();
  await page.setViewportSize({ width: 1000, height: 500 });
  await page.locator('#new-run').click();
  await page.locator('#load-example').click();
  await page.waitForFunction(() => !document.querySelector('#edit-plot').disabled);
  assert.equal(await page.evaluate(() => document.querySelector('#plot-editor').chart.config().plot.minimap.show), true, 'a cropped plot starts with its minimap enabled');
  assert.equal(await page.locator('.clusterMapMinimap').isVisible(), true);
  assert.equal(await page.evaluate(() => document.querySelector('#plot-editor').chart.canUndo()), false, 'automatic navigation defaults do not create an undo step');
  assert.deepEqual(errors, []);
});
