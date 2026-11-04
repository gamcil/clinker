import test from 'node:test';
import assert from 'node:assert/strict';
import { openApp } from './helpers.mjs';

test('Clinker accepts browser File objects from an external opener', { timeout: 30000 }, async t => {
  const { page, url } = await openApp(t);
  await page.goto(url);
  const exampleUrl = new URL('examples/A. versicolor CBS 583.65.gbk', url).href;
  await page.evaluate(async exampleUrl => {
    const response = await fetch(exampleUrl);
    const file = new File([await response.arrayBuffer()], 'handoff.gbk', { type: 'text/plain' });
    window.__handoffResult = null;
    window.addEventListener('message', event => {
      if (event.source !== window.__clinkerPopup || event.origin !== location.origin) return;
      if (event.data?.type === 'clinker:ready') {
        event.source.postMessage({
          type: 'clinker:import-files', requestId: 'browser-test', files: [file],
        }, event.origin);
      } else if (event.data?.type === 'clinker:files-received') {
        window.__handoffResult = event.data;
      }
    });
    const button = document.createElement('button');
    button.id = 'open-integration';
    button.textContent = 'Open clinker';
    button.addEventListener('click', () => {
      const target = new URL(location.href);
      target.search = '?integration=postmessage';
      window.__clinkerPopup = window.open(target.href, 'clinker-integration-test');
    });
    document.body.append(button);
  }, exampleUrl);
  const popupPromise = page.waitForEvent('popup');
  await page.locator('#open-integration').click();
  const clinker = await popupPromise;
  await clinker.waitForLoadState();
  await clinker.waitForFunction(() => {
    const summary = document.querySelector('#file-input-status-summary');
    return summary?.textContent === '1 input file ready';
  });

  assert.equal(await clinker.locator('#file-input-status-list li').count(), 1);
  assert.match(await clinker.locator('#file-input-status-list li').textContent(), /handoff\.gbk/);
  await page.waitForFunction(() => window.__handoffResult?.ok === true);
  assert.equal(await page.evaluate(() => window.__handoffResult.requestId), 'browser-test');
  await clinker.close();
});
