import { createHash } from 'node:crypto';
import { mkdir, readdir, readFile, rename, rm, writeFile } from 'node:fs/promises';
import { resolve, relative } from 'node:path';
import { assets, copy, repoRoot, webRoot } from './assets.mjs';

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

async function files(directory) {
  const result = [];
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    if (entry.isDirectory()) result.push(...await files(path));
    else result.push(path);
  }
  return result.sort();
}

// Every input participates, including vendor chunks, WASM, examples and builder.
const hash = createHash('sha256');
for (const path of await files(staging)) {
  hash.update(relative(staging, path)).update('\0').update(await readFile(path)).update('\0');
}
for (const name of ['package-lock.json', 'scripts/assets.mjs', 'scripts/build.mjs']) {
  hash.update(name).update('\0').update(await readFile(resolve(webRoot, name))).update('\0');
}
const id = hash.digest('hex').slice(0, 16);
for (const path of await files(staging)) {
  const name = relative(staging, path);
  // Upstream modules keep their own content-hashed chunk paths intact.
  if (name.startsWith('vendor/') || name.startsWith('pkg/') || !/\.(js|mjs|html)$/.test(name)) continue;
  let text = await readFile(path, 'utf8');
  text = text.replace(/(["'])((?:\.\/)?(?:vendor\/)?[\w./-]+\.(?:js|mjs|css|wasm))\1/g,
    (_, quote, url) => `${quote}${url}?v=${id}${quote}`);
  await writeFile(path, text);
}
const inventory = {};
for (const path of await files(staging)) {
  inventory[relative(staging, path)] = createHash('sha256').update(await readFile(path)).digest('hex');
}
await writeFile(resolve(staging, 'build.json'), JSON.stringify({ format: 'clinker-web-build', version: 1, id, files: inventory }, null, 2) + '\n');
await rm(resolve(webRoot, 'dist'), { recursive: true, force: true });
await rename(staging, resolve(webRoot, 'dist'));
console.log(`Built web/dist (${id}).`);
