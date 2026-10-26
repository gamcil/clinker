import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { initSync, BrowserFileParser, AlignmentWorker, post_process, open_session } from '../../pkg/clinker_wasm.js';

initSync({ module: readFileSync(new URL('../../pkg/clinker_wasm_bg.wasm', import.meta.url)) });
const fixture = name => readFileSync(new URL(`../../../crates/clinker-core/tests/fixtures/${name}`, import.meta.url));
const settings = { enabled: false, kmerSize: 5, minSharedKmers: 8, identityCutoff: 0.3, genbank: { ranges: [], asSeparateClusters: false, setOrigin: true } };
const genbank = (name, translation, record = 'TEST') => ({ name, bytes: new TextEncoder().encode(`LOCUS       ${record.padEnd(26)}42 bp    DNA     linear   UNA 01-JAN-2000\nFEATURES             Location/Qualifiers\n     CDS             1..42\n                     /locus_tag="gene"\n                     /translation="${translation}"\nORIGIN\n        1 atgatgatgatgatgatgatgatgatgatgatgatgatgatg\n//\n`) });
function parse(files, options = settings) {
  const parser = new BrowserFileParser(options);
  try {
    for (const file of files) parser.parse_file(file);
    return parser.finish(() => {});
  } finally { parser.free(); }
}

test('incremental WASM parsing and post-processing preserve mixed-format input order', () => {
  const parsed = parse([
    { name: 'features.gff3', bytes: fixture('features.gff3') },
    genbank('last.gbk', 'MSTAVKDLGYRQAA'),
    { name: 'features.fna', bytes: fixture('features.fna') },
  ]);
  assert.deepEqual(parsed.layout.clusters.map(c => c.source_file), ['features.gff3', 'last.gbk']);
  const result = post_process(parsed.layout, parsed.proteins, [], [], 0.3, true, () => {});
  assert.deepEqual(result.plotData.clusters.map(c => c.source_file), ['features.gff3', 'last.gbk']);
  const reopened = open_session(JSON.parse(JSON.stringify(result.analysisSession)));
  assert.deepEqual(reopened.similarityMatrix, result.similarityMatrix);
  assert.deepEqual(reopened.analysisSession, result.analysisSession);
});

test('prefilter and alignment agree on case and unsupported amino acids', () => {
  const files = [genbank('one.gbk', 'MOUSTAVKDLGYRQ'), genbank('two.gbk', 'mxxstavkdlgyrq')];
  for (const kmerSize of [5, 6, 7, 8]) {
    const parsed = parse(files, { ...settings, enabled: true, kmerSize, minSharedKmers: kmerSize === 5 ? 8 : 2 });
    assert.deepEqual(parsed.candidatePairs, [{ queryIndex: 0, targetIndex: 1 }]);
  }
  const parsed = parse(files);
  const worker = new AlignmentWorker();
  try {
    const matches = worker.analyse_batch(parsed.proteins, parsed.alignmentBatches[0], [], 0.3, () => {});
    assert.equal(matches.links.length, 1);
    assert.equal(matches.links[0].identity, 1);
  } finally { worker.free(); }
});

test('explicit circular trim retains original coordinates despite automatic rotation', () => {
  const parsed = parse([{ name: 'circular.gbk', bytes: fixture('circular.gbk') }], {
    ...settings, genbank: { ...settings.genbank, ranges: [{ record: 'CIRCULAR', start: 9, end: 18 }] },
  });
  const locus = parsed.layout.clusters[0].loci[0];
  assert.deepEqual([locus.start, locus.end], [9, 18]);
  assert.deepEqual(locus.genes.map(g => [g.label, g.start, g.end]), [['middle', 9, 18]]);
});

test('record splitting options reach the incremental parser', () => {
  const bytes = Buffer.concat([genbank('x', 'MSTAVKDLGYRQAA', 'FIRST').bytes, genbank('x', 'MSTAVKDLGYRQAA', 'SECOND').bytes]);
  const parsed = parse([{ name: 'records.gbk', bytes }], { ...settings, genbank: { ...settings.genbank, asSeparateClusters: true } });
  assert.deepEqual(parsed.layout.clusters.map(c => c.name), ['FIRST', 'SECOND']);
});

test('incremental parsing applies frame, split CDS phase and genetic code to derived proteins', () => {
  const gff = '##gff-version 3\nseq\t.\tCDS\t1\t1\t.\t+\t2\tID=split;transl_table=4\nseq\t.\tCDS\t2\t11\t.\t+\t1\tID=split;transl_table=4\n##FASTA\n>seq\nCCATGTGATAA\n';
  const parsed = parse([
    { name: 'translation.gbk', bytes: fixture('translation.gbk') },
    { name: 'partial.gff3', bytes: new TextEncoder().encode(gff) },
  ]);
  assert.deepEqual(parsed.proteins.map(p => p.translation), ['MW', 'MW', 'MW']);
  const result = post_process(parsed.layout, parsed.proteins, [], [], 0.3, true, () => {});
  assert.deepEqual(result.analysisSession.analysis.clusters.flatMap(c => c.loci.flatMap(l => l.genes.map(g => g.translation))), ['MW', 'MW', 'MW']);
});


test('post-processing retains Ward merges only when similarity ordering runs', () => {
  const parsed = parse([genbank('a.gbk', 'MSTAVKDLGYRQAA'), genbank('b.gbk', 'MSTAVKDLGYRQAA'), genbank('c.gbk', 'MSTAVKDLGYRQAA')]);
  const ref = cluster => ({ cluster, locus: 0, gene: 0 });
  const links = [{ query: ref(0), target: ref(1), identity: 0.95, similarity: 0.95 }, { query: ref(1), target: ref(2), identity: 0.7, similarity: 0.7 }];
  const result = post_process(parsed.layout, parsed.proteins, links, [], 0.3, false, () => {});
  assert.equal(result.clusteringMerges.length, 2);
  assert.equal(result.clusteringHeights.length, 2);
  assert.ok(result.clusteringHeights.every(height => Number.isFinite(height) && height > 0));
  const leaves = [];
  const stack = [4];
  while (stack.length) {
    const node = stack.pop();
    if (node < 3) leaves.push(node);
    else stack.push(...[...result.clusteringMerges[node - 3]].reverse());
  }
  assert.deepEqual(leaves, result.clusterOrder);
  assert.deepEqual(post_process(parsed.layout, parsed.proteins, links, [], 0.3, true, () => {}).clusteringMerges, []);
});


test('alignment workers reuse buffers for successive Cartesian and selected batches', () => {
  const worker = new AlignmentWorker();
  try {
    for (const length of [3, 100, 5]) {
      const protein = { cluster: 0, locus: 0, gene: 0, translation: 'MSTAVKLMDEFGH'.repeat(8).slice(0, length) };
      const target = { ...protein, cluster: 1 };
      const tile = worker.analyse_batch([protein, target], { kind: 'cartesian', query: { start: 0, end: 1 }, target: { start: 1, end: 2 } }, [], 0.3, () => {});
      const pairs = worker.analyse_batch([protein, target], { kind: 'selected', pairs: { start: 0, end: 1 } }, [{ queryIndex: 0, targetIndex: 1 }], 0.3, () => {});
      assert.deepEqual(pairs, tile);
      assert.equal(tile.alignmentCount, 1);
      assert.equal(tile.links.length, 1);
      assert.equal(tile.links[0].identity, 1);
    }
    const parsed = parse([genbank('a.gbk', 'MSTAVKLMDEFGH'), genbank('b.gbk', 'MSTAVKLMDEFGH')]);
    assert.deepEqual(parsed.alignmentBatches, [{ kind: 'cartesian', query: { start: 0, end: 1 }, target: { start: 1, end: 2 } }]);
    assert.deepEqual(parsed.alignmentCounts, { possiblePairs: 1, plannedPairs: 1, batchCount: 1 });
    assert.throws(() => worker.analyse_batch(parsed.proteins, { kind: 'selected', pairs: { start: 0, end: 1 } }, [], 0.3, () => {}), /out of bounds/);
    const filtered = parse([genbank('a.gbk', 'MSTAVKLMDEFGH'), genbank('b.gbk', 'MSTAVKLMDEFGH')], { ...settings, enabled: true });
    assert.deepEqual(filtered.alignmentBatches, [{ kind: 'selected', pairs: { start: 0, end: 1 } }]);
    assert.deepEqual(filtered.alignmentCounts, { possiblePairs: 1, plannedPairs: 1, batchCount: 1 });
    const rejected = parse([genbank('a.gbk', 'MSTAVKLMDEFGH'), genbank('b.gbk', 'QQQQQQQQQQQQQ')], { ...settings, enabled: true });
    assert.deepEqual(rejected.alignmentCounts, { possiblePairs: 1, plannedPairs: 0, batchCount: 0 });
  } finally { worker.free(); }
});
