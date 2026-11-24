import test from 'node:test';
import assert from 'node:assert/strict';
import { openApp } from './helpers.mjs';

test('the quick tour is optional and can be reopened from help', { timeout: 90000 }, async t => {
  const { page, url } = await openApp(t);
  await page.goto(url);
  assert.equal(await page.locator('#kmer-size').isVisible(), false);
  await page.locator('.run-other-settings > summary').click();
  assert.equal(await page.locator('#kmer-size').isVisible(), true);
  await page.locator('#prefilter').uncheck();
  assert.equal(await page.locator('#kmer-size').isVisible(), false);
  await page.locator('#prefilter').check();
  assert.equal(await page.locator('#kmer-size').isVisible(), true);
  await page.locator('.run-other-settings > summary').click();
  await page.locator('#load-example').click();
  await page.locator('#tour-prompt').waitFor({ state: 'visible' });
  const original = await page.evaluate(() => {
    const chart = document.querySelector('#plot-editor').chart;
    chart.patch([{ type: 'groups.update', ids: [chart.data().groups[0].uid], changes: { label: 'Edited group' } }]);
    window.originalTourChart = chart;
    return { state: chart.state(), data: chart.data(), selection: chart.selection(), canUndo: chart.canUndo(), canRedo: chart.canRedo() };
  });
  await page.getByRole('button', { name: 'Quick tour', exact: true }).first().click();
  await page.getByText('New run', { exact: true }).last().waitFor();
  for (const title of ['Save and export', 'Add data', 'Row order', 'Figure editor', 'Edit figure data', 'Adjust appearance', 'Theme', 'Help', 'Navigate', 'Select loci', 'Arrange loci', 'Arrange clusters', 'Tree actions', 'Legend']) {
    await page.locator('.driver-popover-next-btn').click();
    await page.locator('.driver-popover-title').getByText(title, { exact: true }).waitFor();
    if (title === 'Select loci') {
      assert.equal(await page.locator('.tour-locus-target').count(), 1);
      assert.notEqual(await page.locator('.driver-overlay').evaluate(el => getComputedStyle(el).opacity), '0');
      assert.ok(await page.locator('.driver-popover kbd').count() >= 5);
      assert.equal(await page.evaluate(() => {
        const chart = document.querySelector('#plot-editor').chart;
        const ids = chart.visibleData().clusters.flatMap(cluster => cluster.loci.map(locus => locus.uid));
        const bounds = chart.screenBounds({ loci: ids });
        const target = document.querySelector('.tour-locus-target').getBoundingClientRect();
        const viewport = document.querySelector('#plot').getBoundingClientRect();
        return Math.abs(target.top - Math.max(bounds.top, viewport.top)) < 1 &&
          Math.abs(target.bottom - Math.min(bounds.bottom, viewport.bottom)) < 1;
      }), true);
      const { uid, target } = await page.evaluate(() => {
        const chart = document.querySelector('#plot-editor').chart;
        const uid = chart.visibleData().clusters[0].loci[0].uid;
        return { uid, target: chart.screenBounds({ loci: [uid] }) };
      });
      await page.keyboard.down('Shift');
      await page.mouse.click(target.left + target.width * 0.25, target.top + target.height / 2);
      await page.keyboard.up('Shift');
      await page.waitForFunction(id => document.querySelector('#plot-editor').chart.selection().ids.includes(id), uid);
      const toolbar = page.locator('.tour-example-toolbar');
      await toolbar.waitFor({ state: 'visible' });
      await page.waitForFunction(() => document.querySelector('.tour-example-toolbar').classList.contains('driver-allowed-element'));
      await toolbar.getByRole('button', { name: 'Clear', exact: true }).click();
      assert.deepEqual(await page.evaluate(() => document.querySelector('#plot-editor').chart.selection().ids), []);
      const gene = await page.evaluate(() => {
        const chart = document.querySelector('#plot-editor').chart;
        return chart.screenBounds({ genes: [chart.visibleData().clusters[0].loci[0].genes[0].uid] });
      });
      await page.mouse.click(gene.left + gene.width / 2, gene.top + gene.height / 2, { button: 'right' });
      const geneMenu = page.locator('.tour-example .tooltip');
      await geneMenu.getByText('Edit label', { exact: true }).waitFor({ state: 'visible' });
      assert.equal(await geneMenu.evaluate(el => el.classList.contains('driver-allowed-element')), true);
      assert.ok(await geneMenu.evaluate(el => Number(getComputedStyle(el).zIndex)) >
        await page.locator('.driver-overlay').evaluate(el => Number(getComputedStyle(el).zIndex)));
      await geneMenu.locator('input[type="text"]').first().fill('Tour gene');
      await page.evaluate(() => { document.documentElement.dataset.theme = 'dark'; });
      assert.equal(await page.locator('.driver-popover').evaluate(el => getComputedStyle(el).backgroundColor), 'rgb(35, 42, 53)');
      assert.equal(await page.locator('.driver-popover kbd').first().evaluate(el => getComputedStyle(el).color), 'rgb(229, 234, 241)');
      await page.evaluate(() => { document.documentElement.dataset.theme = 'light'; });
    }
    if (title === 'Arrange clusters') {
      await page.locator('.tour-locus-target').waitFor({ state: 'detached' });
      assert.equal(await page.locator('.tour-region-target[data-region="labels"]').count(), 1);
    }
    if (title === 'Tree actions') {
      assert.equal(await page.locator('.tour-region-target[data-region="tree"]').count(), 1);
      assert.ok(await page.evaluate(() => document.querySelector('#plot-editor').chart.data().tree));
      assert.equal(await page.evaluate(() => {
        const chart = document.querySelector('#plot-editor').chart;
        const scale = chart.screenBounds({ treeScaleBar: true });
        const target = document.querySelector('.tour-region-target[data-region="tree"]').getBoundingClientRect();
        const viewport = document.querySelector('#plot').getBoundingClientRect();
        return scale && target.bottom >= Math.min(scale.bottom, viewport.bottom) - 1;
      }), true);
      const node = await page.evaluate(() => document.querySelector('#plot-editor').chart.screenBounds({ treeNodes: ['clustering-5'] }));
      await page.mouse.click(node.left + node.width / 2, node.top + node.height / 2, { button: 'right' });
      const menu = page.locator('.tour-example .tooltip');
      await menu.getByRole('button', { name: 'Rotate branches' }).waitFor({ state: 'visible' });
      assert.ok(await menu.evaluate(el => Number(getComputedStyle(el).zIndex)) >
        await page.locator('.driver-overlay').evaluate(el => Number(getComputedStyle(el).zIndex)));
      assert.equal(await page.locator('.tour-example-chart').evaluate(el => getComputedStyle(el).isolation), 'auto');
    }
    if (title === 'Arrange clusters' || title === 'Tree actions') {
      assert.equal(await page.evaluate(() => {
        const chart = document.querySelector('#plot-editor').chart;
        const ids = chart.visibleData().clusters.flatMap(cluster => cluster.loci.map(locus => locus.uid));
        const bounds = chart.screenBounds({ loci: ids });
        const target = document.querySelector('.tour-region-target').getBoundingClientRect();
        const viewport = document.querySelector('#plot').getBoundingClientRect();
        return target.top <= Math.max(bounds.top, viewport.top) + 1 &&
          target.bottom >= Math.min(bounds.bottom, viewport.bottom) - 1;
      }), true);
    }
    if (title === 'Legend') {
      const target = page.locator('.tour-region-target[data-region="legend"]');
      assert.equal(await target.count(), 1);
      const bounds = await target.boundingBox();
      const prompt = await page.locator('.driver-popover').boundingBox();
      assert.ok(prompt.y + prompt.height <= bounds.y);
      await page.mouse.click(bounds.x + 30, bounds.y + 10, { button: 'right' });
      const menu = page.locator('.tour-example .tooltip');
      await menu.locator('input').first().waitFor({ state: 'visible' });
      const after = await page.locator('.driver-popover').boundingBox();
      assert.equal(after.x, prompt.x);
      assert.equal(after.y, prompt.y);
      assert.notEqual(await page.locator('.driver-overlay').evaluate(el => getComputedStyle(el).opacity), '0');
      assert.equal(await menu.evaluate(el => el.classList.contains('driver-allowed-element')), true);
      assert.ok(await menu.evaluate(el => Number(getComputedStyle(el).zIndex)) >
        await page.locator('.driver-overlay').evaluate(el => Number(getComputedStyle(el).zIndex)));
      assert.equal(await menu.locator('input').first().evaluate(el => {
        const box = el.getBoundingClientRect();
        return document.elementFromPoint(box.left + box.width / 2, box.top + box.height / 2) === el;
      }), true);
      await menu.locator('input').first().fill('Tour group');
    }
    if (title === 'Figure editor') assert.equal(await page.locator('.tour-region-target').count(), 0);
    if (title === 'Navigate') {
      await page.locator('.tour-locus-target').waitFor({ state: 'detached' });
      const locus = await page.evaluate(() => document.querySelector('#plot-editor').chart.visibleData().clusters.find(cluster => cluster.uid === 'cluster-3').loci[0]);
      assert.equal(locus.bio.start, 7571);
    }
    if (title === 'Edit figure data' || title === 'Adjust appearance') {
      assert.equal(await page.locator('#plot-editor').isVisible(), true);
      const active = title === 'Edit figure data' ? 'Data' : 'Appearance';
      assert.equal(await page.locator('#plot-editor').getByRole('tab', { name: active, exact: true }).getAttribute('aria-selected'), 'true');
    }
    if (title === 'Theme') assert.equal(await page.locator('#plot-editor').isHidden(), true);
  }
  await page.locator('.driver-popover-close-btn').click();
  assert.equal(await page.locator('#tour-prompt').isHidden(), true);
  assert.deepEqual(await page.evaluate(() => {
    const chart = document.querySelector('#plot-editor').chart;
    return { state: chart.state(), data: chart.data(), selection: chart.selection(), canUndo: chart.canUndo(), canRedo: chart.canRedo() };
  }), original);
  assert.equal(await page.evaluate(() => document.querySelector('#plot-editor').chart === window.originalTourChart), true);
  assert.equal(original.canUndo, true);
  assert.equal(await page.locator('.tour-example').count(), 0);
  await page.locator('#plot-help').click();
  await page.locator('#help-tour').click();
  await page.locator('.driver-popover-title').getByText('New run', { exact: true }).waitFor();
  for (const title of ['Save and export', 'Add data', 'Row order', 'Figure editor', 'Edit figure data', 'Adjust appearance', 'Theme', 'Help', 'Navigate', 'Select loci']) {
    await page.locator('.driver-popover-next-btn').click();
    await page.locator('.driver-popover-title').getByText(title, { exact: true }).waitFor();
  }
  await page.locator('.driver-popover-title').getByText('Select loci', { exact: true }).waitFor();
  assert.equal(await page.locator('.tour-locus-target').count(), 1);
  await page.locator('.driver-popover-close-btn').click();
  assert.equal(await page.locator('.tour-locus-target').count(), 0);
  assert.deepEqual(await page.evaluate(() => document.querySelector('#plot-editor').chart.state().camera), original.state.camera);
  await page.setViewportSize({ width: 390, height: 844 });
  await page.evaluate(() => {
    const chart = document.querySelector('#plot-editor').chart;
    chart.config({ plot: { renderer: 'canvas', minimap: { show: true } } });
    chart.selection({ kind: 'loci', ids: [chart.visibleData().clusters[0].loci[0].uid] });
  });
  await page.waitForFunction(() => {
    const toolbar = document.querySelector('.cm-plot-selection-toolbar');
    const minimap = document.querySelector('canvas.clusterMapMinimap');
    if (!toolbar || toolbar.hidden || !minimap) return false;
    const a = toolbar.getBoundingClientRect(), b = minimap.getBoundingClientRect();
    return a.right <= b.left || a.left >= b.right || a.bottom <= b.top || a.top >= b.bottom;
  });
  await page.locator('#plot-help').click();
  await page.locator('#help-tour').click();
  await page.locator('.driver-popover-title').getByText('New run', { exact: true }).waitFor();
  await page.evaluate(() => { document.querySelector('#plot-editor').chart.screenBounds = () => null; });
  const fallbackCamera = await page.evaluate(() => document.querySelector('#plot-editor').chart.state().camera);
  for (const title of ['Save and export', 'Add data', 'Row order', 'Figure editor', 'Edit figure data', 'Adjust appearance', 'Theme', 'Help', 'Navigate', 'Select loci']) {
    await page.locator('.driver-popover-next-btn').click();
    await page.locator('.driver-popover-title').getByText(title, { exact: true }).waitFor();
  }
  assert.equal(await page.locator('.tour-locus-target').count(), 0);
  assert.deepEqual(await page.evaluate(() => document.querySelector('#plot-editor').chart.state().camera), fallbackCamera);
  await page.locator('.driver-popover-close-btn').click();
  assert.equal(await page.locator('.tour-locus-target').count(), 0);
});
