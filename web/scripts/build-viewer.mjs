import { mkdir, readFile, rm, writeFile } from 'node:fs/promises';
import { resolve } from 'node:path';
import { fileURLToPath } from 'node:url';
import { assets, copy, repoRoot, webRoot } from './assets.mjs';
import { finalizeBuild } from './build-common.mjs';

export async function buildViewer() {
  const staging = resolve(webRoot, '.build/viewer-site');
  await rm(staging, { recursive: true, force: true });
  await mkdir(staging, { recursive: true });
  // Share the shell; the native entry has no run form or analysis entry script.
  const shell = (await readFile(resolve(webRoot, 'src/index.html'), 'utf8'))
    .replace(/\s*<button id="new-run"[^>]*>.*?<\/button>/s, '')
    .replace(/\s*<!-- run-panel:start -->[\s\S]*?<!-- run-panel:end -->/, '')
    .replace('src="app.js"', 'src="local.js"');
  if (shell.includes('id="run-panel"')) throw new Error('Shared shell is missing run-panel markers.');
  await writeFile(resolve(staging, 'index.html'), shell);
  for (const name of ['local.js', 'style.css', 'viewer', 'tree-import.mjs', 'tree-input.mjs', 'annotation-import.mjs', 'workflow.mjs']) {
    await copy(resolve(webRoot, 'src', name), resolve(staging, name));
  }
  for (const asset of assets) await copy(resolve(webRoot, asset.source), resolve(staging, asset.output));
  await copy(resolve(repoRoot, 'LICENSE'), resolve(staging, 'LICENSE'));
  await finalizeBuild(staging, resolve(webRoot, '.build/viewer'), 'viewer', ['scripts/build-viewer.mjs']);
}

if (process.argv[1] === fileURLToPath(import.meta.url)) await buildViewer();
