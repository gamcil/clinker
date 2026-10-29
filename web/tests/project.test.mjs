import { parseTwoColumnCsv } from "../src/delimited.mjs";
import { parseAnnotationTable } from "../src/annotation-import.mjs";
import test from "node:test";
import assert from "node:assert/strict";
import { similarityCsv, parseProject, displayedClusterOrder } from "../src/viewer/project.mjs";

test("CSV follows edited names, order and deleted clusters while retaining analysis indices", () => {
  const pair = { queryCoverage: 0.8, targetCoverage: 0.6, similarity: 0.8 };
  const matrix = Array.from({ length: 3 }, () => [pair, pair, pair]);
  const project = { data: { clusters: ['A, one', 'B', 'C'].map((name, index) => ({ uid: `cluster-${index}`, name })) }, state: { clusterOrder: [] } };
  assert.deepEqual(displayedClusterOrder(project), ['cluster-0', 'cluster-1', 'cluster-2']);
  project.state.clusterOrder = ['cluster-2', 'cluster-0', 'cluster-1'];
  assert.equal(similarityCsv(matrix, project).split('\n')[1], '2,"A, one",0,3,B,1,0.8,0.6,0.8,0.19999999999999996');
  project.data.clusters.splice(1, 1);
  project.data.clusters[0].name = 'Renamed';
  project.state.clusterOrder.pop();
  const csv = similarityCsv(matrix, project).trim().split('\n');
  assert.equal(csv.length, 2);
  assert.equal(csv[1], '2,Renamed,0,1,C,2,0.8,0.6,0.8,0.19999999999999996');
  project.data.clusters[0].uid = 'missing';
  assert.throws(() => similarityCsv(matrix, project), /cluster IDs/);
});

function project() {
  const locus = { uid: 'locus-0-0', start: 0, end: 30, genes: [
    { uid: 'gene-0-0-0', start: 0, end: 9 },
    { uid: 'gene-0-0-1', start: 10, end: 20 },
  ] };
  return { format: 'clinker-project', version: 2, config: {},
    analysis: { analysis: { clusters: [{ loci: [{ genes: [{}, {}] }] }] } },
    data: { clusters: [{ uid: 'cluster-0', loci: [locus] }], links: [], groups: [] },
    state: { version: 2, clusterOrder: [] },
  };
}

test("project validation permits edits and deletions but rejects broken IDs, coordinates and layout", () => {
  const saved = project();
  saved.data.clusters[0].name = 'Edited';
  saved.data.clusters[0].loci[0].genes.pop();
  saved.data.links.push({ uid: 'link-0', query: { uid: 'gene-0-0-0' }, target: { uid: 'gene-0-0-1' } });
  assert.deepEqual(parseProject(JSON.stringify(saved)), saved);
  const cases = [
    [value => value.data.clusters.push(value.data.clusters[0]), /unique uid/],
    [value => { value.data.clusters[0].uid = 'cluster-9'; }, /analysis clusters/],
    [value => { value.data.clusters[0].loci[0].uid = 'locus-1-0'; }, /analysis loci/],
    [value => { value.data.clusters[0].loci[0].genes[0].uid = 'gene-0-0-9'; }, /analysis genes/],
    [value => { value.data.clusters[0].loci[0].genes[0].end = -1; }, /coordinates/],
    [value => { value.data.links = [{ uid: 'link-0', query: { uid: 'missing' }, target: { uid: 'missing' } }]; }, /known genes/],
    [value => { value.state = []; }, /layout state/],
    [value => { value.state.version = 9; }, /layout state/],
    [value => { value.state.clusterOrder = ['cluster-0', 'cluster-0']; }, /cluster order/],
    [value => { value.state.clusterOrder = ['cluster-9']; }, /cluster order/],
  ];
  for (const [mutate, error] of cases) {
    const invalid = project();
    mutate(invalid);
    assert.throws(() => parseProject(JSON.stringify(invalid)), error);
  }
  assert.throws(() => parseProject('null'), /not a clinker project/);
});

test("legacy layout states upgrade gene keys while preserving other edits", () => {
  const saved = project();
  saved.state = { version: 1, clusterOrder: [], camera: { x: 10, y: 20, k: 2 },
    loci: [['locus-0-0', { start: 2, end: 25, flipped: true }]],
    genes: [['locus-0-0:gene-0-0-0', { start: 3, end: 9, strand: -1 }]] };
  const opened = parseProject(JSON.stringify(saved));
  assert.deepEqual(opened.state, { ...saved.state, version: 2,
    genes: [['gene-0-0-0', { start: 3, end: 9, strand: -1 }]] });
  assert.equal(saved.state.version, 1);
  assert.deepEqual(parseProject(JSON.stringify(opened)), opened);
});

test("CSV imports share quoted-field parsing and retain their format validation", () => {
  const rows = '\uFEFF"gene,1","enzyme ""A"""\r\n"gene2","line1\nline2"\r\n';
  const expected = [["gene,1", 'enzyme "A"'], ["gene2", "line1\nline2"]];
  assert.deepEqual(parseTwoColumnCsv(rows, "Gene functions"), expected);
  assert.deepEqual(parseAnnotationTable(`source_file,label\r\n${rows.replace(/^\uFEFF/, "")}`, ",").records, expected);
  assert.throws(() => parseTwoColumnCsv('gene,"bad"text', "Gene functions"), /after a quoted field/);
  assert.throws(() => parseTwoColumnCsv("gene,", "Gene functions"), /two non-empty columns/);
  assert.throws(() => parseAnnotationTable("name,NAME\na,b", ","), /unique/);
});
