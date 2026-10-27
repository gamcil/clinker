import { createHash } from 'node:crypto';
import { readdir, readFile, rename, rm, writeFile } from 'node:fs/promises';
import { resolve, relative } from 'node:path';
import { webRoot } from './assets.mjs';

export async function files(directory) {
  const result = [];
  for (const entry of await readdir(directory, { withFileTypes: true })) {
    const path = resolve(directory, entry.name);
    if (entry.isDirectory()) result.push(...await files(path));
    else result.push(path);
  }
  return result.sort();
}

export async function finalizeBuild(staging, destination, target, builders) {
  const nameOf = path => relative(staging, path).replaceAll('\\', '/');
  const hash = createHash('sha256');
  for (const path of await files(staging)) {
    hash.update(nameOf(path)).update('\0').update(await readFile(path)).update('\0');
  }
  for (const name of ['package-lock.json', 'scripts/assets.mjs', 'scripts/build-common.mjs', ...builders]) {
    hash.update(name).update('\0').update(await readFile(resolve(webRoot, name))).update('\0');
  }
  const id = hash.digest('hex').slice(0, 16);
  for (const path of await files(staging)) {
    const name = nameOf(path);
    if (name.startsWith('vendor/') || name.startsWith('pkg/') || !/\.(js|mjs|html)$/.test(name)) continue;
    const text = (await readFile(path, 'utf8')).replace(
      /(["'])((?:\.\/)?(?:vendor\/)?[\w./-]+\.(?:js|mjs|css|wasm))\1/g,
      (_, quote, url) => `${quote}${url}?v=${id}${quote}`,
    );
    await writeFile(path, text);
  }
  const inventory = {};
  for (const path of await files(staging)) {
    inventory[nameOf(path)] = createHash('sha256').update(await readFile(path)).digest('hex');
  }
  await writeFile(resolve(staging, 'build.json'), JSON.stringify({ format: 'clinker-web-build', version: 1, target, id, files: inventory }, null, 2) + '\n');
  await rm(destination, { recursive: true, force: true });
  await rename(staging, destination);
  console.log(`Built ${target} (${id}).`);
}
