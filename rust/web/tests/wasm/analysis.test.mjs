import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { initSync, BrowserFileParser, analyse_tile, post_process, open_session } from '../../pkg/clinker_wasm.js';

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
  const reopened = open_session(JSON.parse(JSON.stringify(result.analysisSession)), true, () => {});
  assert.deepEqual(reopened.plotData, result.plotData);
});

test('prefilter and alignment agree on case and unsupported amino acids', () => {
  const files = [genbank('one.gbk', 'MOUSTAVKDLGYRQ'), genbank('two.gbk', 'mxxstavkdlgyrq')];
  for (const kmerSize of [5, 6, 7, 8]) {
    const parsed = parse(files, { ...settings, enabled: true, kmerSize, minSharedKmers: kmerSize === 5 ? 8 : 2 });
    assert.deepEqual(parsed.candidatePairs, [{ queryIndex: 0, targetIndex: 1 }]);
  }
  const parsed = parse(files);
  const matches = analyse_tile([parsed.proteins[0]], [parsed.proteins[1]], 0.3, () => {});
  assert.equal(matches.links.length, 1);
  assert.equal(matches.links[0].identity, 1);
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
