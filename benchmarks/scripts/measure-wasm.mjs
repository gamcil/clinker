#!/usr/bin/env node
import { chromium } from 'playwright';
import { readdir, writeFile } from 'node:fs/promises';
import { spawn } from 'node:child_process';
import { createServer } from 'node:net';
import { basename, resolve } from 'node:path';
import process from 'node:process';

function argumentsFrom(argv) {
  const values = {};
  for (let index = 0; index < argv.length; index += 2) {
    values[argv[index].replace(/^--/, '').replaceAll('-', '_')] = argv[index + 1];
  }
  return values;
}

async function availablePort() {
  const server = createServer();
  await new Promise((resolveListen, reject) => {
    server.once('error', reject);
    server.listen(0, '127.0.0.1', resolveListen);
  });
  const { port } = server.address();
  await new Promise(resolveClose => server.close(resolveClose));
  return port;
}

async function waitForServer(url, serverProcess) {
  for (let attempt = 0; attempt < 100; attempt++) {
    if (serverProcess.exitCode !== null) throw new Error('Static web server exited before it became ready.');
    try {
      const response = await fetch(url);
      if (response.ok) return;
    } catch {}
    await new Promise(resolveWait => setTimeout(resolveWait, 100));
  }
  throw new Error('Timed out waiting for the local web server.');
}

const args = argumentsFrom(process.argv.slice(2));
const mode = args.mode || 'wasm';
if (!['wasm', 'wasm-kmer'].includes(mode)) throw new Error(`Invalid mode: ${mode}`);
const fixture = resolve(args.fixture);
const output = resolve(args.output);
const stem = `${output}`;
const logPath = `${stem}.log`;
const jsonPath = `${stem}.json`;
const figurePath = `${stem}.figure.json`;
const inputs = (await readdir(fixture))
  .filter(name => /\.(gb|gbk|genbank|gbf|gbff)$/i.test(name))
  .sort()
  .map(name => resolve(fixture, name));
if (!inputs.length) throw new Error(`No GenBank files found in ${fixture}`);

const port = await availablePort();
const server = spawn('python3', [
  '-m', 'http.server', String(port), '--bind', '127.0.0.1', '--directory', '/opt/clinker-web',
], { stdio: 'ignore' });
let browser;
const log = [];

try {
  const appUrl = `http://127.0.0.1:${port}/`;
  await waitForServer(appUrl, server);
  browser = await chromium.launch({ headless: true });
  const page = await browser.newPage({ viewport: { width: 1280, height: 800 } });
  await page.addInitScript(threads => {
    Object.defineProperty(Navigator.prototype, 'hardwareConcurrency', {
      configurable: true,
      get: () => threads,
    });
  }, Number(args.threads));
  page.on('console', message => {
    if (message.type() === 'error') log.push(`console: ${message.text()}`);
  });
  page.on('pageerror', error => log.push(`pageerror: ${error.stack || error}`));
  await page.route('**/benchmark-input', route => route.fulfill({
    contentType: 'text/html', body: '<!doctype html><input id="files" type="file" multiple>',
  }));
  await page.goto(`${appUrl}benchmark-input`);
  await page.locator('#files').setInputFiles(inputs);
  await page.evaluate(async () => {
    window.benchmarkModules = await Promise.all([
      import('/analysis.mjs'), import('/alignment.mjs'), import('/workflow.mjs'),
    ]);
  });
  const result = await page.evaluate(async ({ mode, identity, threads, kmerSize, shared }) => {
    const [{ parseInputFiles, postProcess }, { analyseTilesInWorkerPool }, { Workflow }] = window.benchmarkModules;
    const workflow = new Workflow();
    const settings = {
      enabled: mode === 'wasm-kmer', noAlign: false, identityCutoff: identity,
      kmerSize, minSharedKmers: shared,
      genbank: { ranges: [], asSeparateClusters: false, setOrigin: true },
    };
    const started = performance.now();
    const files = [];
    for (const file of document.querySelector('#files').files) {
      files.push({ name: file.name, bytes: new Uint8Array(await workflow.readFile(file)) });
    }
    const parsed = await parseInputFiles(files, settings, {
      onProgress() {}, onPrefilterStart() {}, onPrefilterProgress() {},
    }, workflow);
    const aligned = await analyseTilesInWorkerPool(parsed, identity, parsed.wasmModule, () => {}, workflow, threads);
    const prepared = await postProcess({
      layout: parsed.layout, proteins: parsed.proteins, links: aligned.links,
      groupDefinitions: [], settings, useFileOrder: false, wasmModule: parsed.wasmModule,
    }, () => {}, workflow);
    const figure = JSON.stringify(prepared.plotData);
    const jsonBytes = new TextEncoder().encode(figure).byteLength;
    const elapsedSeconds = (performance.now() - started) / 1000;
    workflow.dispose();
    return { figure, jsonBytes, elapsedSeconds, retainedLinks: prepared.plotData.links.length,
      groups: prepared.plotData.groups.length, alignmentCount: aligned.alignmentCount };
  }, { mode, identity: Number(args.identity), threads: Number(args.threads),
       kmerSize: Number(args.kmer_size), shared: Number(args.min_shared_kmers) });
  if (log.some(line => line.startsWith('pageerror:'))) throw new Error('The web app reported a browser error.');

  const record = {
    dataset: basename(fixture),
    role: args.role,
    clusters: Number(args.clusters),
    proteins: Number(args.proteins),
    mode,
    repeat: Number(args.repeat),
    identity: Number(args.identity),
    threads: Number(args.threads),
    kmer_size: mode === 'wasm-kmer' ? Number(args.kmer_size) : null,
    min_shared_kmers: mode === 'wasm-kmer' ? Number(args.min_shared_kmers) : null,
    browser_version: browser.version(),
    elapsed_seconds: result.elapsedSeconds,
    timing_boundary: 'figure-json-v1',
    retained_links: result.retainedLinks,
    groups: result.groups,
    json_bytes: result.jsonBytes,
    alignments_performed: result.alignmentCount,
  };
  await writeFile(figurePath, result.figure);
  await writeFile(jsonPath, `${JSON.stringify(record, null, 2)}\n`);
  log.push(JSON.stringify(record));
  process.stdout.write(`${JSON.stringify(record)}\n`);
} catch (error) {
  log.push(error.stack || String(error));
  throw error;
} finally {
  await browser?.close();
  server.kill('SIGTERM');
  await writeFile(logPath, `${log.join('\n')}\n`);
}
