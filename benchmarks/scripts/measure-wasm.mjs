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
if (args.timing_boundary !== 'retained-links-v1') throw new Error('Invalid timing boundary.');
const fixture = resolve(args.fixture);
const output = resolve(args.output);
const stem = `${output}`;
const logPath = `${stem}.log`;
const jsonPath = `${stem}.json`;
const linksPath = `${stem}.links`;
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
  await page.goto(appUrl);
  await page.locator('#identity').fill(String(Number(args.identity) * 100));
  if (mode === 'wasm-kmer') {
    await page.locator('#prefilter').check();
    await page.locator('#kmer-size').fill(args.kmer_size);
    await page.locator('#min-shared-kmers').fill(args.min_shared_kmers);
  } else {
    await page.locator('#prefilter').uncheck();
  }
  await page.locator('#files').setInputFiles(inputs);
  await page.waitForFunction(() => !document.querySelector('#analyse').disabled);

  await page.evaluate(() => {
    const step = document.querySelector('[data-run-step="align"]');
    const result = step.querySelector('.step-result');
    const textContent = Object.getOwnPropertyDescriptor(Node.prototype, 'textContent');
    // Capture the completion synchronously, before figure preparation begins.
    Object.defineProperty(result, 'textContent', {
      configurable: true,
      get() { return textContent.get.call(this); },
      set(value) {
        textContent.set.call(this, value);
        if (step.dataset.state === 'complete' && !window.__benchmarkAlignment) {
          window.__benchmarkAlignment = {
            elapsedSeconds: (performance.now() - window.__benchmarkStartedAt) / 1000,
            summary: value,
          };
        }
      },
    });
    window.__benchmarkStartedAt = performance.now();
    document.querySelector('#analyse').click();
  });
  await page.waitForFunction(() => {
    const editor = document.querySelector('#plot-editor');
    const panel = document.querySelector('#run-panel');
    const data = editor?.chart?.data?.();
    return panel?.hidden && Array.isArray(data?.links);
  }, null, { timeout: 30 * 60 * 1000 });
  const result = await page.evaluate(() => {
    const links = document.querySelector('#plot-editor').chart.data().links;
    const alignment = window.__benchmarkAlignment;
    if (!alignment) throw new Error('Alignment completion was not recorded.');
    return {
      elapsedSeconds: alignment.elapsedSeconds,
      links,
      alignmentSummary: alignment.summary.trim(),
    };
  });
  if (log.some(line => line.startsWith('pageerror:'))) throw new Error('The web app reported a browser error.');

  const alignmentCount = result.alignmentSummary.match(/^Ran ([\d,]+) alignments/);
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
    timing_boundary: args.timing_boundary,
    retained_links: result.links.length,
    alignments_performed: alignmentCount ? Number(alignmentCount[1].replaceAll(',', '')) : null,
    alignment_summary: result.alignmentSummary,
  };
  await writeFile(linksPath, `${JSON.stringify(result.links)}\n`);
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
