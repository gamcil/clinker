import test from 'node:test';
import assert from 'node:assert/strict';
import { createServer } from 'node:http';
import { mkdir, readFile } from 'node:fs/promises';
import { resolve, extname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

const root = resolve(process.env.CLINKER_WEB_ROOT || fileURLToPath(new URL('../../', import.meta.url)));
const browserType = { chromium, firefox, webkit }[process.env.CLINKER_TEST_BROWSER || 'chromium'];

test('task menus support adding data, saving, editing and keyboard navigation at narrow widths', { timeout: 90000 }, async t => {
  const server = createServer(async (request, response) => {
    const path = resolve(root, '.' + decodeURIComponent(new URL(request.url, 'http://localhost').pathname).replace(/\/$/, '/index.html'));
    if (!path.startsWith(root + '/')) { response.writeHead(403).end(); return; }
    try {
      response.setHeader('Content-Type', { '.html': 'text/html', '.js': 'text/javascript', '.mjs': 'text/javascript', '.css': 'text/css', '.wasm': 'application/wasm' }[extname(path)] || 'application/octet-stream');
      response.end(await readFile(path));
    } catch { response.writeHead(404).end(); }
  });
  await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
  t.after(() => { server.closeAllConnections(); return new Promise(resolve => server.close(resolve)); });
  const browser = await browserType.launch();
  t.after(() => browser.close());
  const page = await browser.newPage({ viewport: { width: 1400, height: 900 } });
  const errors = [];
  page.on('pageerror', error => errors.push(error.message));
  await page.goto(`http://127.0.0.1:${server.address().port}/`);
  await page.waitForFunction(() => Boolean(document.querySelector('#run-tree-file')));
  const openMenu = id => page.locator(`${id} > summary`).click();
  const capture = async name => {
    if (!process.env.CLINKER_NAV_SCREENSHOTS) return;
    await mkdir(process.env.CLINKER_NAV_SCREENSHOTS, { recursive: true });
    await page.screenshot({ path: resolve(process.env.CLINKER_NAV_SCREENSHOTS, `${name}.png`) });
  };

  await page.locator('#file-menu > summary').focus();
  await page.keyboard.press('ArrowDown');
  assert.equal(await page.locator('#load-clinker-file').evaluate(button => button === document.activeElement), true);
  assert.equal(await page.locator('#download-project').isEnabled(), false);
  await page.keyboard.press('Escape');
  assert.equal(await page.locator('#file-menu').evaluate(menu => menu.open), false);
  await page.waitForFunction(() => document.activeElement === document.querySelector('#file-menu > summary'));
  assert.equal(await page.locator('#run-panel').isVisible(), true);

  await page.locator('#load-example').click();
  await page.waitForFunction(() => !document.querySelector('#edit-plot').disabled);
  const menuGaps = await page.locator('#file-menu > summary, #add-data-menu > summary, #plot-order > summary').evaluateAll(items => {
    const boxes = items.map(item => item.getBoundingClientRect());
    return [boxes[1].left - boxes[0].right, boxes[2].left - boxes[1].right];
  });
  assert.ok(Math.abs(menuGaps[0] - menuGaps[1]) < 1, 'the workflow menus have consistent spacing');
  await capture('desktop-light');
  await openMenu('#add-data-menu');
  await page.locator('#import-annotations').click();
  assert.equal(await page.locator('#annotation-dialog').isVisible(), true, errors.join('\n'));
  await page.locator('#annotation-dialog').getByRole('button', { name: 'Cancel', exact: true }).click();
  await page.waitForFunction(() => document.activeElement === document.querySelector('#add-data-menu > summary'));

  await openMenu('#add-data-menu');
  await page.locator('#manage-tree').click();
  await page.locator('#tree-file').setInputFiles(resolve(root, 'examples/bundled-clusters.nwk'));
  await page.waitForFunction(() => !document.querySelector('#apply-tree').disabled);
  await page.locator('#apply-tree').click();
  assert.equal(await page.locator('#plot-order-label').textContent(), 'Tree');
  await openMenu('#plot-order');
  await page.locator('[data-order="clustering"]').click();
  assert.equal(await page.locator('#plot-order-label').textContent(), 'Similarity');

  await openMenu('#plot-order');
  await page.locator('#tree-appearance').click();
  assert.equal(await page.locator('#edit-plot').textContent(), 'Close editor');
  assert.equal(await page.locator('#plot-editor').isVisible(), true);
  const appearance = page.locator('#plot-editor clinker-appearance-panel');
  const treeField = appearance.locator('[data-config="tree.show"]');
  assert.equal(await page.locator('#plot-editor').getByRole('tab', { name: 'Appearance', exact: true }).getAttribute('aria-selected'), 'true');
  const waitForTreeFocus = () => page.waitForFunction(() => {
    const panel = document.querySelector('#plot-editor').shadowRoot.querySelector('clinker-appearance-panel');
    return panel.shadowRoot.activeElement?.dataset.config === 'tree.show';
  });
  await waitForTreeFocus();
  await capture('tree-settings');
  const filter = appearance.getByRole('searchbox', { name: 'Find an appearance setting', exact: true });
  await filter.fill('legend');
  assert.equal(await treeField.isVisible(), false);
  await page.locator('#edit-plot').click();
  await openMenu('#plot-order');
  await page.locator('#tree-appearance').click();
  assert.equal(await filter.inputValue(), '', 'tree settings are revealed even after filtering other appearance settings');
  await waitForTreeFocus();
  assert.equal(await treeField.evaluate(input => {
    const field = input.getBoundingClientRect();
    const viewport = input.closest('.settings').getBoundingClientRect();
    return field.top >= viewport.top && field.bottom <= viewport.bottom;
  }), true, 'the tree controls are scrolled into view');
  await page.locator('#edit-plot').click();
  assert.equal(await page.locator('#plot-editor').isVisible(), false);

  await openMenu('#file-menu');
  const saving = page.waitForEvent('download');
  await page.locator('#download-project').click();
  const project = JSON.parse(await readFile(await (await saving).path(), 'utf8'));
  assert.equal(project.data.clusters.length, 5);
  assert.ok(project.analysis);
  await page.locator('#new-run').click();
  assert.equal(await page.locator('#run-panel').isVisible(), true);
  await page.locator('#file-menu > summary').focus();
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('End');
  assert.equal(await page.locator('#download-matrix').evaluate(button => button === document.activeElement), true);
  await page.keyboard.press('Tab');
  await page.waitForFunction(() => !document.querySelector('#file-menu').open);
  await page.locator('#file-menu > summary').focus();
  await page.keyboard.press('ArrowDown');
  await page.keyboard.press('Escape');
  assert.equal(await page.locator('#run-panel').isVisible(), true, 'closing a menu does not dismiss comparison setup');
  await page.locator('#cancel-run').click();

  await openMenu('#appearance-menu');
  await page.locator('[data-theme-mode="dark"]').click();
  await capture('desktop-dark');
  await openMenu('#file-menu');
  assert.equal(await page.locator('#download-project').evaluate(button => getComputedStyle(button).borderStyle), 'none');
  await capture('desktop-file-dark');
  await page.keyboard.press('Escape');

  for (const width of [900, 375, 320]) {
    await page.setViewportSize({ width, height: 850 });
    const boxes = await page.locator('.app-bar > h1, #new-run, #file-menu > summary, .app-actions > .menu > summary, #edit-plot, #appearance-label, #plot-help')
      .evaluateAll(items => items.map(item => {
        const { x, y, width, height } = item.getBoundingClientRect();
        return { name: item.textContent, x, y, width, height };
      }));
    for (const box of boxes) {
      assert.ok(box.x >= 0 && box.x + box.width <= width, `${box.name} fits at ${width}px`);
      for (const other of boxes.filter(other => other !== box)) {
        const overlap = Math.min(box.x + box.width, other.x + other.width) - Math.max(box.x, other.x) > 1
          && Math.min(box.y + box.height, other.y + other.height) - Math.max(box.y, other.y) > 1;
        assert.equal(overlap, false, `${box.name} and ${other.name} do not overlap at ${width}px`);
      }
    }
    if (width <= 375) {
      const actions = await page.locator('#add-data-menu > summary, #edit-plot').evaluateAll(items => items.map(item => item.getBoundingClientRect().y));
      assert.ok(Math.abs(actions[0] - actions[1]) < 2, 'plot actions fit on one row below the project controls');
    }
    for (const id of ['#file-menu', '#add-data-menu', '#plot-order', '#appearance-menu']) {
      await openMenu(id);
      const box = await page.locator(`${id} .menu-popover`).boundingBox();
      assert.ok(box.x >= 0 && box.x + box.width <= width, `${id} stays inside ${width}px viewport`);
      await page.keyboard.press('Escape');
    }
    await capture(`narrow-dark-${width}`);
  }
  await page.setViewportSize({ width: 1400, height: 900 });
  await openMenu('#appearance-menu');
  await page.locator('[data-theme-mode="auto"]').click();
  await page.emulateMedia({ colorScheme: 'light' });
  await page.waitForFunction(() => document.documentElement.dataset.theme === 'light');
  assert.deepEqual(errors, []);
});
