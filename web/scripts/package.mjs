import { createHash } from 'node:crypto';
import { lstat, mkdir, mkdtemp, readdir, readFile, rename, rm } from 'node:fs/promises';
import { dirname, resolve, sep } from 'node:path';
import { copy, repoRoot, webRoot } from './assets.mjs';

if (process.argv.length !== 3) throw new Error('Usage: npm --prefix web run package -- /path/to/site');
const destination = resolve(process.argv[2]);
const source = resolve(webRoot, 'dist');
const within = (path, parent) => path === parent || path.startsWith(parent + sep);
const protectedRoots = ['web', 'crates', 'examples', 'docs', 'bench', '.git'].map(name => resolve(repoRoot, name));
if (within(repoRoot, destination) || protectedRoots.some(path => within(destination, path))) {
  throw new Error('Choose an output directory outside repository sources.');
}
const manifest = JSON.parse(await readFile(resolve(source, 'build.json'), 'utf8'));
if (manifest.format !== 'clinker-web-build' || manifest.version !== 1) throw new Error('Build the web app first.');
// Package exactly the completed build. Extra files in dist cannot leak into it.
for (const [name, digest] of Object.entries(manifest.files)) {
  const path = resolve(source, name);
  if (!within(path, source)) throw new Error(`Invalid build path: ${name}`);
  const actual = createHash('sha256').update(await readFile(path)).digest('hex');
  if (actual !== digest) throw new Error(`Build asset changed: ${name}. Rebuild before packaging.`);
}
const existing = await lstat(destination).catch(error => {
  if (error.code === 'ENOENT') return null;
  throw error;
});
if (existing) {
  if (!existing.isDirectory() || existing.isSymbolicLink()) throw new Error('Output must be a directory, not a symlink.');
  if ((await readdir(destination)).length) {
    const previous = await readFile(resolve(destination, 'build.json'), 'utf8').catch(() => '{}');
    if (JSON.parse(previous).format !== 'clinker-web-build') {
      throw new Error('Output is not a previous clinker build. Choose a new or empty directory.');
    }
  }
}
await mkdir(dirname(destination), { recursive: true });
const staging = await mkdtemp(resolve(dirname(destination), '.clinker-package-'));
try {
  for (const name of [...Object.keys(manifest.files), 'build.json']) {
    await copy(resolve(source, name), resolve(staging, name));
  }
  await rm(destination, { recursive: true, force: true });
  await rename(staging, destination);
} finally {
  await rm(staging, { recursive: true, force: true });
}
console.log(`Packaged ${manifest.id} to ${destination}.`);
