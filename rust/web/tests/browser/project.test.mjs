import test from 'node:test';
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium } from 'playwright';

const root = resolve(process.env.CLINKER_WEB_ROOT || fileURLToPath(new URL('../../', import.meta.url)));

test('analysis → edit → project → reopen → CSV and SVG preserves the displayed figure', { timeout: 120000 }, async t => {
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
  const browser = await chromium.launch({ headless: true });
  t.after(() => browser.close());
  const page = await browser.newPage();
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  await page.locator('#files').setInputFiles([
    resolve(root, 'examples/A. versicolor CBS 583.65.gbk'),
    resolve(root, 'examples/A. burnettii MST-FP2249.gbk'),
  ]);
  await page.waitForFunction(() => !document.querySelector('#analyse').disabled);
  await page.locator('#analyse').click();
  await page.waitForFunction(() => document.querySelector('#run-panel').hidden, {}, { timeout: 60000 });
  const expected = await page.evaluate(() => {
    const chart = document.querySelector('#plot-editor').chart;
    const gene = chart.data().clusters[0].loci[0].genes[0];
    chart.patch([{ type: 'genes.update', ids: [gene.uid], changes: { label: 'Audit edited gene' } }]);
    chart.config({ plot: { renderer: 'svg' }, gene: { label: { show: true } } });
    const state = chart.state();
    state.clusterOrder = chart.data().clusters.map(c => c.uid).reverse();
    chart.state(state);
    return { gene: gene.uid, state: chart.state() };
  });
  async function download(id) {
    await page.getByText('Export', { exact: true }).click();
    const downloading = page.waitForEvent('download');
    await page.locator(id).click();
    return readFile(await (await downloading).path(), 'utf8');
  }
  function checkCsv(csv) {
    for (const line of csv.trim().split('\n').slice(1)) {
      const row = line.split(',');
      assert.equal(Number(row[0]), expected.state.clusterOrder.indexOf(`cluster-${row[2]}`) + 1);
      assert.equal(Number(row[3]), expected.state.clusterOrder.indexOf(`cluster-${row[5]}`) + 1);
    }
  }
  checkCsv(await download('#download-matrix'));
  const project = await download('#download-project');
  await page.reload();
  await page.locator('#clinker-file').setInputFiles({ name: 'edited.json', mimeType: 'application/json', buffer: Buffer.from(project) });
  await page.waitForFunction(() => document.querySelector('#plot-status').textContent.includes('Loaded clinker project'), {}, { timeout: 60000 });
  const reopened = await page.evaluate(() => document.querySelector('#plot-editor').chart.project());
  assert.deepEqual(reopened.state.clusterOrder, expected.state.clusterOrder);
  assert.equal(reopened.data.clusters.flatMap(c => c.loci.flatMap(l => l.genes)).find(g => g.uid === expected.gene).label, 'Audit edited gene');
  checkCsv(await download('#download-matrix'));
  assert.match(await download('#download-svg'), /Audit edited gene/);
  assert.deepEqual(errors, []);
});
