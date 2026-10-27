import { createServer } from 'node:http';
import { readFile } from 'node:fs/promises';
import { resolve, extname } from 'node:path';
import { fileURLToPath } from 'node:url';
import { chromium, firefox, webkit } from 'playwright';

export const root = resolve(process.env.CLINKER_WEB_ROOT || fileURLToPath(new URL('../../dist/', import.meta.url)));

export async function openApp(t, pageOptions) {
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
  const browserName = process.env.CLINKER_TEST_BROWSER || 'chromium';
  const browserType = { chromium, firefox, webkit }[browserName];
  if (!browserType) throw new Error(`Unknown test browser: ${browserName}`);
  const browser = await browserType.launch();
  t.after(() => browser.close());
  return { page: await browser.newPage(pageOptions), url: `http://127.0.0.1:${server.address().port}/` };
}
