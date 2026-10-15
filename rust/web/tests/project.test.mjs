import test from "node:test";
import assert from "node:assert/strict";
import { similarityCsv, parseProject, displayedClusterOrder } from "../project.mjs";

test("CSV positions follow displayed stable IDs, including after a saved reorder", () => {
  const pair = { queryCoverage: 0.8, targetCoverage: 0.6, similarity: 0.8 };
  const result = { clusterNames: ['A, one', 'B', 'C'], similarityMatrix: Array.from({ length: 3 }, () => [pair, pair, pair]) };
  const project = { data: { clusters: [0, 1, 2].map(index => ({ uid: `cluster-${index}` })) }, state: { clusterOrder: [] } };
  assert.deepEqual(displayedClusterOrder(project), ['cluster-0', 'cluster-1', 'cluster-2']);
  project.state.clusterOrder = ['cluster-2', 'cluster-0', 'cluster-1'];
  const csv = similarityCsv(result, displayedClusterOrder(JSON.parse(JSON.stringify(project))));
  assert.equal(csv.split('\n')[1], '2,"A, one",0,3,B,1,0.8,0.6,0.8,0.19999999999999996');
  assert.throws(() => similarityCsv(result, ['cluster-0']), /cluster IDs/);
});

test("project loading rejects broken cross references before replacing the chart", () => {
  const saved = { format: 'clinker-project', version: 2, config: {}, analysis: {}, data: { clusters: [], links: [{ query: { uid: 'missing' }, target: { uid: 'missing' } }], groups: [] } };
  assert.throws(() => parseProject(JSON.stringify(saved)), /known genes/);
  assert.throws(() => parseProject('null'), /not a clinker project/);
});
