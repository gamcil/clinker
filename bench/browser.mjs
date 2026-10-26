import { chromium } from '../web/node_modules/playwright/index.mjs';
import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname, sep } from 'node:path';
import { fileURLToPath } from 'node:url';
import readline from 'node:readline';

const root = resolve(fileURLToPath(new URL('..', import.meta.url)));
let inputs = [];
const server = createServer(async (req, res) => {
  try {
    const url = new URL(req.url, 'http://localhost');
    if (url.pathname === '/') {
      res.setHeader('Content-Type', 'text/html');
      res.end('<!doctype html><script type="module" src="/bench/wasm.mjs"></script>');
      return;
    }
    let path;
    if (/^\/input\/\d+$/.test(url.pathname)) path = inputs[Number(url.pathname.split('/').pop())];
    else {
      path = resolve(root, '.' + decodeURIComponent(url.pathname));
      if (!path.startsWith(root + sep) || !['/bench/', '/web/'].some(prefix => url.pathname.startsWith(prefix))) {
        res.writeHead(403).end(); return;
      }
    }
    if (!path) { res.writeHead(404).end(); return; }
    res.setHeader('Content-Type', { '.mjs': 'text/javascript', '.js': 'text/javascript', '.wasm': 'application/wasm' }[extname(path)] || 'application/octet-stream');
    res.end(await readFile(path));
  } catch { res.writeHead(404).end(); }
});
await new Promise(resolve => server.listen(0, '127.0.0.1', resolve));
const browser = await chromium.launch({ headless: !process.argv.includes('--headed') });
process.once('SIGTERM', () => { browser.close().finally(() => process.exit(0)); });
console.log(JSON.stringify({ ready: true, browser: browser.version(), node: process.version }));
try {
  for await (const line of readline.createInterface({ input: process.stdin })) {
    const options = JSON.parse(line);
    inputs = options.files;
    const page = await browser.newPage();
    try {
      await page.goto(`http://127.0.0.1:${server.address().port}/`);
      await page.waitForFunction(() => typeof window.runBenchmark === 'function');
      console.log(JSON.stringify(await page.evaluate(options => window.runBenchmark(options), options)));
    } catch (error) { console.log(JSON.stringify({ error: String(error) })); }
    finally { await page.close(); }
  }
} finally {
  await browser.close();
  server.closeAllConnections();
  await new Promise(resolve => server.close(resolve));
}
