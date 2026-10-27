import { mkdir, readdir, rm } from 'node:fs/promises';
import { resolve } from 'node:path';
import { assets, copy, repoRoot, webRoot } from './assets.mjs';
import { finalizeBuild } from './build-common.mjs';
import { buildViewer } from './build-viewer.mjs';

const staging = resolve(webRoot, '.build/site');
await rm(staging, { recursive: true, force: true });
await mkdir(staging, { recursive: true });
await copy(resolve(webRoot, 'src'), staging);
await copy(resolve(webRoot, '.build/pkg'), resolve(staging, 'pkg'));
for (const asset of assets) await copy(resolve(webRoot, asset.source), resolve(staging, asset.output));

// Copy only the public fixtures, never local benchmark/download directories.
const examples = resolve(repoRoot, 'examples');
for (const entry of await readdir(examples)) {
  if (entry.endsWith('.gbk') || entry === 'bundled-clusters.nwk') {
    await copy(resolve(examples, entry), resolve(staging, 'examples', entry));
  }
}
for (const entry of await readdir(resolve(examples, 'gff3'))) {
  if (/\.(gff3|fna)$/.test(entry)) {
    await copy(resolve(examples, 'gff3', entry), resolve(staging, 'examples/gff3', entry));
  }
}
await copy(resolve(repoRoot, 'LICENSE'), resolve(staging, 'LICENSE'));
await copy(resolve(webRoot, 'DEPLOYMENT.md'), resolve(staging, 'README.md'));

await finalizeBuild(staging, resolve(webRoot, 'dist'), 'app', ['scripts/build.mjs']);
await buildViewer();
