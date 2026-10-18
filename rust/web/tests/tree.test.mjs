import test from 'node:test';
import assert from 'node:assert/strict';
import { readFileSync } from 'node:fs';
import { parseNewick, matchTreeLeaves, resolveTree, validateTree, inputTreeClusters, bindInputTree } from '../tree.mjs';
import { parseProject, displayedClusterOrder } from '../project.mjs';

const clusters = [
  { uid: 'cluster-0', name: 'First cluster', source_file: 'sample_A.gbk', loci: [{ name: 'ABC.1' }] },
  { uid: 'cluster-1', name: 'Second cluster', source_file: 'sample_B.gff3', loci: [{ name: 'DEF.2' }, { name: 'scaffold_2' }] },
  { uid: 'cluster-2', name: 'Third cluster', source_file: 'sample_C.gb', loci: [{ name: 'GHI.1' }] },
];

test('the bundled illustrative tree is fully bifurcating', () => {
  const tree = parseNewick(readFileSync(new URL('../../../examples/bundled-clusters.nwk', import.meta.url), 'utf8'));
  assert.equal(tree.nodes.filter(node => !node.children).length, 5);
  assert.ok(tree.nodes.filter(node => node.children).every(node => node.children.length === 2));
});

test('Newick preserves quoted, numeric and duplicate labels, child order and branch metadata', () => {
  const tree = parseNewick(" [&R] (('A, B':0,'O''Brien':1e-3)95:2E+1,123:0.5,123)named_root; ", { supportScale: 'percent' });
  const root = tree.nodes[0];
  const internal = tree.nodes[1];
  assert.equal(tree.rooted, true);
  assert.deepEqual(root.children, ['tree-1', 'tree-4', 'tree-5']);
  assert.equal(root.label, 'named_root');
  assert.equal(internal.label, '95');
  assert.deepEqual(internal.support, { value: 95, scale: 'percent' });
  assert.equal(internal.branchLength, 20);
  assert.deepEqual(tree.nodes.filter(node => !node.children).map(node => node.label), ['A, B', "O'Brien", '123', '123']);
  assert.equal(tree.nodes[2].branchLength, 0);
  assert.equal(tree.nodes[3].branchLength, 0.001);
  assert.equal(tree.nodes[5].branchLength, undefined);
  assert.equal(tree.nodes[4].support, undefined);
  assert.equal(parseNewick('[&U] (A,B);').rooted, false);
  assert.equal(parseNewick('((A,B)95,C);', { internalLabels: 'names' }).nodes[1].support, undefined);
  const deep = parseNewick('('.repeat(1500) + 'A' + ')'.repeat(1500) + ';');
  assert.equal(deep.nodes.length, 1501);
});

test('Newick rejects invalid lengths, unsupported annotations, malformed syntax and implicit root metadata loss', () => {
  for (const length of ['-1', 'NaN', 'Infinity', '1e999', "'0.1'"]) {
    assert.throws(() => parseNewick(`(A:${length},B);`), /finite, non-negative/);
  }
  for (const input of ['', '(A,B)', '(A,B;', 'A,B;', '(A,,B);', '(A B);', '(A,B));', '(A,B);(C,D);', "('A,B);", '(A:1:2,B);']) {
    assert.throws(() => parseNewick(input), Error, input);
  }
  assert.throws(() => parseNewick('(A[&&NHX:S=x],B);'), /annotation/);
  assert.throws(() => parseNewick('(A[comment],B);'), /annotation/);
  assert.throws(() => parseNewick('((A,B)101,C);', { supportScale: 'percent' }), /support/);
  assert.throws(() => parseNewick('((A,B)95,C);', { supportScale: 'proportion' }), /support/);
  assert.throws(() => parseNewick('((A,B)1e999,C);'), /support/);
  assert.throws(() => parseNewick('(A,B):0;'), /display root/);
  assert.throws(() => parseNewick('(A,B)100;'), /display root/);
  const tree = parseNewick('(A,B)100:0;', { ignoreRootMetadata: true });
  assert.equal(tree.nodes[0].label, '100');
  assert.equal(tree.nodes[0].branchLength, undefined);
  assert.equal(tree.nodes[0].support, undefined);
});

test('tree matching uses exact metadata, detects collisions and supports explicit one-to-one corrections', () => {
  const tree = parseNewick('(sample_C,(ABC.1,sample_B.gff3));');
  const resolved = resolveTree(tree, clusters);
  assert.deepEqual(validateTree(resolved, new Set(clusters.map(c => c.uid))), ['cluster-2', 'cluster-0', 'cluster-1']);
  assert.equal(tree.nodes.find(node => node.label === 'ABC.1').cluster, undefined, 'resolution does not mutate parsed input');
  for (const name of ['abc.1', 'ABC', 'sample_a', 'First_cluster']) {
    assert.equal(matchTreeLeaves(parseNewick(`('${name}',sample_B,sample_C);`), clusters).ready, false, name);
  }
  const extra = matchTreeLeaves(parseNewick('(sample_A,sample_B,sample_C,extra);'), clusters);
  assert.match(extra.rows.at(-1).error, /No match/);
  const incomplete = matchTreeLeaves(parseNewick('(sample_A,sample_B);'), clusters);
  assert.deepEqual(incomplete.missing.map(c => c.uid), ['cluster-2']);
  const repeated = parseNewick('(sample_A,sample_A,sample_C);');
  assert.match(matchTreeLeaves(repeated, clusters).rows[0].error, /more than once/);
  const assignments = new Map([[repeated.nodes[2].uid, 'cluster-1']]);
  assert.equal(matchTreeLeaves(repeated, clusters, assignments).ready, true);
  assert.doesNotThrow(() => resolveTree(repeated, clusters, assignments));
  const ambiguous = clusters.map(c => ({ ...c, name: 'Shared name' }));
  const result = matchTreeLeaves(parseNewick("('Shared name',sample_B,sample_C);"), ambiguous);
  assert.match(result.rows[0].error, /Ambiguous/);
  assert.throws(() => resolveTree(parseNewick('(missing,sample_B,sample_C);'), clusters), /Match every/);
});

test('project trees retain metadata and authoritative child order, rejecting invalid graphs before mounting', () => {
  const tree = resolveTree(parseNewick('((sample_C:0,sample_A:1)0.95:2,sample_B:3);'), clusters);
  const data = { clusters: clusters.map(c => ({ ...c, loci: [] })), links: [], groups: [], tree };
  const order = ['cluster-2', 'cluster-0', 'cluster-1'];
  const saved = { format: 'clinker-project', version: 2, config: {},
    analysis: { analysis: { clusters: clusters.map(() => ({ loci: [] })) } },
    data, state: { version: 1, clusterOrder: order } };
  assert.deepEqual(parseProject(JSON.stringify(saved)), saved);
  assert.deepEqual(displayedClusterOrder({ data, state: { clusterOrder: [] } }), order);
  const cases = [
    value => { value.tree.nodes[0].children.push('missing'); },
    value => { value.tree.nodes[1].children.push(value.tree.root); },
    value => { value.tree.nodes.push({ uid: 'unreachable', cluster: 'cluster-0' }); },
    value => { value.tree.nodes[2].cluster = 'cluster-1'; },
    value => { value.tree.nodes[0].branchLength = 0; },
    value => { value.tree.nodes[1].support.scale = 'unknown'; },
    value => { value.tree.nodes[1].children = []; },
    value => { value.tree.nodes[1].children = null; },
    value => { value.tree.nodes[1].cluster = 'cluster-0'; },
  ];
  for (const mutate of cases) {
    const invalid = structuredClone(saved);
    mutate(invalid.data);
    assert.throws(() => parseProject(JSON.stringify(invalid)), /tree/i);
  }
  saved.state.clusterOrder = [...order].reverse();
  assert.throws(() => parseProject(JSON.stringify(saved)), /tree child order/);
  saved.state.clusterOrder = [];
  saved.inactiveTree = saved.data.tree;
  assert.throws(() => parseProject(JSON.stringify(saved)), /active and inactive/);
  delete saved.data.tree;
  assert.deepEqual(parseProject(JSON.stringify(saved)), saved);
  saved.inactiveTree.nodes[0].children.push('missing');
  assert.throws(() => parseProject(JSON.stringify(saved)), /Invalid tree/);
});


test('input tree matches survive file reordering and bind to final cluster IDs', () => {
  const files = clusters.map(cluster => ({ file: { name: cluster.source_file }, records: cluster.loci }));
  const inputs = inputTreeClusters(files);
  const tree = resolveTree(parseNewick('(sample_A,sample_B,sample_C);'), inputs);
  const reordered = [...clusters].reverse().map((cluster, i) => ({ ...cluster, uid: `final-${i}` }));
  const bound = bindInputTree(tree, reordered);
  assert.deepEqual(validateTree(bound, new Set(reordered.map(c => c.uid))), ['final-2', 'final-1', 'final-0']);
  assert.throws(() => bindInputTree(tree, reordered.slice(1)), /inputs changed/);
});

test('input matching splits GenBank records only and retains sequence and organism aliases', () => {
  const records = [{ name: 'one', organism: 'Organism A' }, { name: 'two', organism: 'Organism B' }];
  const files = [{ file: { name: 'multi.gbk' }, records }, { file: { name: 'multi.gff3' }, records }];
  assert.equal(inputTreeClusters(files).length, 2);
  const split = inputTreeClusters(files, true);
  assert.equal(split.length, 3);
  assert.equal(split[0].name, 'Organism A');
  assert.equal(matchTreeLeaves(parseNewick("('Organism A','Organism B',multi.gff3);"), split).ready, true);
  assert.throws(() => inputTreeClusters([{ file: { name: 'duplicate.gbk' }, records: [records[0], records[0]] }], true), /distinct record/);
});
