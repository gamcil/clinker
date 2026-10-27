// Maintainer/CI check. Cargo consumers need no Node or frontend tooling.
import { spawn, spawnSync } from 'node:child_process';
import { createHash } from 'node:crypto';
import { mkdtemp, readFile, writeFile, rm, realpath } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import { resolve, join } from 'node:path';
import { fileURLToPath } from 'node:url';

const root = fileURLToPath(new URL('../', import.meta.url));
const options = process.argv.slice(2);
if (options.some(arg => !['--offline', '--allow-dirty'].includes(arg))) {
  throw new Error('Usage: node scripts/check-packages.mjs [--offline] [--allow-dirty]');
}
const offline = options.includes('--offline') ? ['--offline'] : [];
const target = resolve(root, 'target/package-check');
function run(command, args, cwd = root, capture = false) {
  const result = spawnSync(command, args, { cwd, encoding: 'utf8', stdio: capture ? 'pipe' : 'inherit' });
  if (result.error) throw result.error;
  if (result.status !== 0) throw new Error(`${command} failed (${result.status})${capture ? `: ${result.stderr}` : ''}`);
  return result.stdout;
}
const metadata = JSON.parse(run('cargo', ['metadata', '--locked', '--no-deps', '--format-version', '1'], root, true));
const cli = metadata.packages.find(pkg => pkg.name === 'clinker-cli');
const core = metadata.packages.find(pkg => pkg.name === 'clinker-core');
const binary = cli.targets.find(target => target.kind.includes('bin')).name;
if (process.env.CLINKER_BINARY && process.env.CLINKER_BINARY !== binary) {
  throw new Error('Workflow CLINKER_BINARY does not match the CLI manifest.');
}
await rm(resolve(target, 'package'), { recursive: true, force: true });
// Verification is performed below on the extracted pair. This also works before
// the core crate is published and avoids registry verification of that peer.
run('cargo', ['package', '--locked', '--no-verify', '-p', core.name, '-p', cli.name, '--target-dir', target, ...options]);

const temporary = await realpath(await mkdtemp(join(tmpdir(), 'clinker-packages-')));
try {
  for (const pkg of [core, cli]) {
    const archive = resolve(target, 'package', `${pkg.name}-${pkg.version}.crate`);
    run('tar', ['-xzf', archive, '-C', temporary]);
    const directory = join(temporary, `${pkg.name}-${pkg.version}`);
    for (const name of ['README.md', 'LICENSE', 'Cargo.lock']) await readFile(join(directory, name));
  }
  const coreDirectory = join(temporary, `${core.name}-${core.version}`);
  const cliDirectory = join(temporary, `${cli.name}-${cli.version}`);
  // Published manifests use a registry dependency. Before the first publication,
  // substitute only the peer archive, never a checkout path. Preserve every
  // third-party locked version while changing core's source to this local path.
  const packagedLock = await readFile(join(cliDirectory, 'Cargo.lock'), 'utf8');
  const coreHash = createHash('sha256').update(await readFile(resolve(target, 'package', `${core.name}-${core.version}.crate`))).digest('hex');
  const coreEntry = packagedLock.split('[[package]]').find(block => /^name = "clinker-core"$/m.test(block));
  if (!coreEntry?.includes(`checksum = "${coreHash}"`)) throw new Error('CLI lockfile does not identify the packaged core archive.');
  const lock = packagedLock.split('[[package]]').map(block => {
    return /^name = "clinker-core"$/m.test(block)
      ? block.replace(/^source = .*\n/gm, '').replace(/^checksum = .*\n/gm, '') : block;
  }).join('[[package]]');
  await writeFile(join(temporary, 'Cargo.lock'), lock);
  await writeFile(join(cliDirectory, 'Cargo.lock'), lock);
  await writeFile(join(temporary, 'Cargo.toml'), `[workspace]\nmembers = [${JSON.stringify(coreDirectory)}, ${JSON.stringify(cliDirectory)}]\nresolver = "2"\n\n[patch.crates-io]\nclinker-core = { path = ${JSON.stringify(coreDirectory)} }\n`);
  // Extracted archives have stable timestamps; rebuild same-version project sources.
  run('cargo', ['clean', '--locked', '-p', core.name, '-p', cli.name, '--target-dir', target, ...offline], temporary);
  run('cargo', ['build', '--locked', '--workspace', '--all-targets', '--target-dir', target, ...offline], temporary);
  const installation = join(temporary, 'installation');
  run('cargo', ['install', '--locked', '--debug', '--path', cliDirectory, '--root', installation,
    '--target-dir', target, '--config', `patch.crates-io.clinker-core.path=${JSON.stringify(coreDirectory)}`, ...offline], temporary);
  const executable = join(installation, 'bin', binary + (process.platform === 'win32' ? '.exe' : ''));
  run(executable, ['--version'], temporary);
  const project = join(temporary, 'project.json');
  const portable = join(temporary, 'figure.html');
  run(executable, [join(coreDirectory, 'tests/fixtures/examples/A. versicolor CBS 583.65.gbk'),
    '--no-align', '--project', project, '--plot', portable], temporary);
  if (!(await readFile(portable, 'utf8')).startsWith('<!doctype html>')) throw new Error('Portable HTML was not produced.');
  const manifest = JSON.parse(await readFile(join(cliDirectory, 'assets/viewer/build.json'), 'utf8'));
  // The installed executable must continue to work with both source trees gone.
  await rm(coreDirectory, { recursive: true });
  await rm(cliDirectory, { recursive: true });
  const child = spawn(executable, ['--open-project', project, '--no-browser'], { cwd: temporary });
  try {
    const url = await new Promise((resolve, reject) => {
      let output = '', errors = '';
      const timer = setTimeout(() => reject(new Error(`Viewer did not start: ${errors}`)), 20000);
      child.stderr.on('data', data => { errors += data; });
      child.stdout.on('data', data => {
        output += data;
        const match = /Viewer: (http:\S+)/.exec(output);
        if (match) { clearTimeout(timer); resolve(match[1]); }
      });
      child.once('error', error => { clearTimeout(timer); reject(error); });
      child.once('exit', code => { clearTimeout(timer); reject(new Error(`Viewer exited ${code}: ${errors}`)); });
    });
    const saved = await (await fetch(new URL('api/project', url))).json();
    if (saved.format !== 'clinker-project' || saved.data.clusters.length !== 1) throw new Error('Native project did not load.');
    for (const [name, expected] of Object.entries(manifest.files)) {
      const response = await fetch(new URL(name, url));
      if (!response.ok) throw new Error(`Missing embedded asset: ${name}`);
      const hash = createHash('sha256').update(Buffer.from(await response.arrayBuffer())).digest('hex');
      if (hash !== expected) throw new Error(`Embedded asset differs: ${name}`);
    }
    console.log('Extracted source packages build, install, export HTML/projects and serve every embedded asset.');
  } finally {
    child.kill();
    if (child.exitCode === null && child.signalCode === null) await new Promise(resolve => child.once('exit', resolve));
  }
} finally {
  await rm(temporary, { recursive: true, force: true });
}
