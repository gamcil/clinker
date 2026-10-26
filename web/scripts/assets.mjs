import { cp, mkdir, readFile } from 'node:fs/promises';
import { dirname, resolve } from 'node:path';
import { fileURLToPath } from 'node:url';

export const webRoot = fileURLToPath(new URL('../', import.meta.url));
export const repoRoot = resolve(webRoot, '..');

// One inventory for browser dependencies and the subset embedded by the CLI.
export const assets = [
  { source: 'node_modules/d3/dist/d3.min.js', output: 'vendor/d3/d3.min.js', embedded: 'd3.min.js' },
  { source: 'node_modules/d3/LICENSE', output: 'LICENSE.d3', embedded: 'LICENSE.d3' },
  { source: 'node_modules/clustermap.js/dist/clustermap.min.js', output: 'vendor/clustermap/clustermap.min.js', embedded: 'clustermap.min.js' },
  { source: 'node_modules/clustermap.js/LICENSE', output: 'LICENSE.clustermap', embedded: 'LICENSE.clustermap' },
  { source: 'node_modules/clustermap.js/dist/editor.mjs', output: 'vendor/clustermap/editor.mjs' },
  { source: 'node_modules/clustermap.js/dist/annotationOverlay.mjs', output: 'vendor/clustermap/annotationOverlay.mjs' },
  { source: 'node_modules/clustermap.js/dist/editor', output: 'vendor/clustermap/editor' },
];

export async function copy(source, destination) {
  await mkdir(dirname(destination), { recursive: true });
  await cp(source, destination, { recursive: true });
}

export async function syncEmbedded({ check = false } = {}) {
  const outdated = [];
  for (const asset of assets.filter(asset => asset.embedded)) {
    const source = resolve(webRoot, asset.source);
    const destination = resolve(repoRoot, 'crates/clinker-cli/assets', asset.embedded);
    if (check) {
      const expected = await readFile(source);
      const actual = await readFile(destination).catch(error => {
        if (error.code === 'ENOENT') return Buffer.alloc(0);
        throw error;
      });
      if (!actual.equals(expected)) outdated.push(asset.embedded);
    } else {
      await copy(source, destination);
    }
  }
  if (outdated.length) {
    throw new Error(`Embedded assets differ: ${outdated.join(', ')}. Run npm --prefix web run sync:assets and commit the result.`);
  }
}
