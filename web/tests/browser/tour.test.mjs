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
  await page.getByRole('button', { name: 'Quick tour', exact: true }).first().click();
  await page.getByText('Navigate', { exact: true }).last().waitFor();
  const position = await page.locator('.driver-popover').boundingBox();
  const viewport = page.viewportSize();
  assert.ok(Math.abs(position.y + position.height / 2 - viewport.height / 2) < 5);
  for (const title of ['Select loci', 'New run', 'Save and export', 'Add data', 'Row order', 'Figure editor', 'Edit figure data', 'Adjust appearance', 'Theme', 'Help']) {
    await page.locator('.driver-popover-next-btn').click();
    await page.locator('.driver-popover-title').getByText(title, { exact: true }).waitFor();
    if (title === 'Select loci') {
      assert.ok(await page.locator('.driver-popover kbd').count() >= 5);
      await page.evaluate(() => { document.documentElement.dataset.theme = 'dark'; });
      assert.equal(await page.locator('.driver-popover').evaluate(el => getComputedStyle(el).backgroundColor), 'rgb(35, 42, 53)');
      assert.equal(await page.locator('.driver-popover kbd').first().evaluate(el => getComputedStyle(el).color), 'rgb(229, 234, 241)');
      await page.evaluate(() => { document.documentElement.dataset.theme = 'light'; });
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
  await page.locator('#plot-help').click();
  await page.locator('#help-tour').click();
  await page.locator('.driver-popover-title').getByText('Navigate', { exact: true }).waitFor();
});
