const numeric = /^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$/;

// A deliberately limited Newick reader: unsupported annotations fail explicitly.
// Iterative construction also handles deep, unbalanced trees without recursion.
export function parseNewick(text, { internalLabels = 'support', supportScale = 'unspecified', ignoreRootMetadata = false } = {}) {
  const tokens = [];
  const fail = (message, position = text.length) => { throw new Error(`${message} (character ${position + 1}).`); };
  for (let i = 0; i < text.length;) {
    const position = i;
    const ch = text[i++];
    if (/\s/.test(ch)) continue;
    if ('(),:;'.includes(ch)) { tokens.push({ kind: ch, position }); continue; }
    if (ch === '[') {
      const end = text.indexOf(']', i);
      if (end < 0) fail('Unclosed annotation', position);
      const flag = text.slice(i, end).trim();
      if (tokens.length || !['&R', '&U'].includes(flag)) fail('Only a leading [&R] or [&U] annotation is supported', position);
      tokens.push({ kind: 'rooted', value: flag === '&R', position });
      i = end + 1;
      continue;
    }
    let value = '';
    if (ch === "'") {
      let closed = false;
      while (i < text.length) {
        const next = text[i++];
        if (next !== "'") value += next;
        else if (text[i] === "'") { value += "'"; i++; }
        else { closed = true; break; }
      }
      if (!closed) fail('Unclosed quoted label', position);
    } else {
      if (ch === ']' || ch === '"') fail('Use single quotes around labels containing spaces or punctuation', position);
      value = ch;
      while (i < text.length && !/[\s(),:;\[\]'"\]]/.test(text[i])) value += text[i++];
    }
    tokens.push({ kind: 'label', value, position, quoted: ch === "'" });
  }
  const tree = { root: '', rooted: null, nodes: [] };
  if (tokens[0]?.kind === 'rooted') tree.rooted = tokens.shift().value;
  const parents = [];
  let current;
  let expected = true;
  let ended = false;
  const add = node => {
    node.uid = `tree-${tree.nodes.length}`;
    if (parents.length) parents.at(-1).children.push(node.uid);
    else if (tree.root) fail('Expected one tree');
    else tree.root = node.uid;
    tree.nodes.push(node);
    current = node;
  };
  for (let i = 0; i < tokens.length; i++) {
    const token = tokens[i];
    if (expected) {
      if (token.kind === '(') { add({ children: [] }); parents.push(current); continue; }
      if (token.kind !== 'label' || !token.value) fail('Expected a leaf label or subtree', token.position);
      add({ label: token.value });
      expected = false;
    } else if (token.kind === ':') {
      if (current.branchLength !== undefined) fail('Duplicate branch length', token.position);
      const length = tokens[++i];
      if (length?.kind !== 'label' || length.quoted || !numeric.test(length.value)
          || !Number.isFinite(Number(length.value)) || Number(length.value) < 0) {
        fail('Branch lengths must be finite, non-negative numbers', length?.position ?? token.position);
      }
      current.branchLength = Number(length.value);
    } else if (token.kind === ',') {
      if (!parents.length) fail('A comma must be inside a subtree', token.position);
      expected = true;
    } else if (token.kind === ')') {
      if (!parents.length) fail('Unexpected closing parenthesis', token.position);
      current = parents.pop();
      if (tokens[i + 1]?.kind === 'label') current.label = tokens[++i].value;
    } else if (token.kind === ';') {
      if (parents.length) fail('Unclosed subtree', token.position);
      if (i !== tokens.length - 1) fail('Only one Newick tree may be imported', token.position);
      ended = true;
    } else fail('Expected a comma, branch length, closing parenthesis or semicolon', token.position);
  }
  if (!ended || expected || !tree.root) fail('A complete Newick tree ending with a semicolon is required');
  if (!['names', 'support'].includes(internalLabels)) fail('Unknown internal-label interpretation');
  if (!['unspecified', 'percent', 'proportion'].includes(supportScale)) fail('Unknown support scale');
  for (const node of tree.nodes) {
    if (node.children && internalLabels === 'support' && numeric.test(node.label || '')) {
      node.support = { value: Number(node.label), scale: supportScale };
    }
    if (node.uid === tree.root && (node.branchLength !== undefined || node.support !== undefined)) {
      if (!ignoreRootMetadata) throw Object.assign(new Error('The display root has a branch length or support value. Select “Discard root branch length and support” to continue.'), { code: 'ROOT_BRANCH_METADATA' });
      delete node.branchLength;
      delete node.support;
    } else if (node.support) checkBranch(node);
  }
  return tree;
}

function checkBranch(node) {
  if (node.branchLength !== undefined && (!Number.isFinite(node.branchLength) || node.branchLength < 0)) {
    throw new Error('Tree branch lengths must be finite and non-negative.');
  }
  if (node.support !== undefined) {
    const { value, scale = 'unspecified', type } = node.support || {};
    if (!Number.isFinite(value) || !['unspecified', 'percent', 'proportion'].includes(scale)
        || (type !== undefined && typeof type !== 'string')
        || (scale === 'percent' && (value < 0 || value > 100))
        || (scale === 'proportion' && (value < 0 || value > 1))) throw new Error('Invalid tree support value or scale.');
  }
}

// Validate persisted/imported trees before replacing a working chart; return tip order.
export function validateTree(tree, clusterIds) {
  const fail = message => { throw new Error(`Invalid tree: ${message}`); };
  if (!tree || !Array.isArray(tree.nodes) || (tree.rooted != null && typeof tree.rooted !== 'boolean')) fail('expected nodes and optional biological rootedness.');
  const nodes = new Map();
  for (const node of tree.nodes) {
    if (!node || typeof node.uid !== 'string' || !node.uid.trim() || nodes.has(node.uid)) fail('node IDs must be unique non-empty strings.');
    nodes.set(node.uid, node);
  }
  if (!nodes.has(tree.root)) fail('missing root.');
  const visited = new Set();
  const assigned = new Set();
  const order = [];
  const stack = [tree.root];
  while (stack.length) {
    const uid = stack.pop();
    const node = nodes.get(uid);
    if (!node || visited.has(uid)) fail('missing child, cycle or shared child.');
    visited.add(uid);
    if (node.label !== undefined && typeof node.label !== 'string') fail('labels must be text.');
    checkBranch(node);
    if (uid === tree.root && (node.branchLength !== undefined || node.support !== undefined)) fail('root cannot have incoming branch metadata.');
    if (node.children !== undefined) {
      if (!Array.isArray(node.children) || !node.children.length || node.cluster !== undefined) fail('internal nodes need children and no cluster.');
      for (let i = node.children.length - 1; i >= 0; i--) stack.push(node.children[i]);
    } else {
      if (!clusterIds.has(node.cluster) || assigned.has(node.cluster)) fail('each leaf must reference a distinct plotted cluster.');
      assigned.add(node.cluster);
      order.push(node.cluster);
    }
  }
  if (visited.size !== nodes.size || assigned.size !== clusterIds.size) fail('all nodes and plotted clusters must occur exactly once.');
  return order;
}

export function matchTreeLeaves(tree, clusters, assignments = new Map()) {
  const aliases = new Map();
  for (const cluster of clusters) {
    const names = [cluster.uid, cluster.name, ...cluster.loci.map(locus => locus.name)];
    if (cluster.source_file) {
      const basename = cluster.source_file.split(/[\\/]/).at(-1);
      names.push(cluster.source_file, basename, basename.replace(/\.(gb|gbk|genbank|gbf|gbff|gff|gff3)$/i, ''));
    }
    for (const name of new Set(names.filter(name => typeof name === 'string' && name.length))) {
      if (!aliases.has(name)) aliases.set(name, []);
      aliases.get(name).push(cluster.uid);
    }
  }
  const rows = tree.nodes.filter(node => !node.children).map(node => {
    const candidates = aliases.get(node.label) || [];
    return { uid: node.uid, label: node.label, candidates,
      cluster: assignments.has(node.uid) ? assignments.get(node.uid) : candidates.length === 1 ? candidates[0] : '' };
  });
  const counts = new Map();
  for (const row of rows) if (row.cluster) counts.set(row.cluster, (counts.get(row.cluster) || 0) + 1);
  const ids = new Set(clusters.map(cluster => cluster.uid));
  for (const row of rows) {
    row.error = !row.cluster ? row.candidates.length > 1 ? 'Ambiguous match' : 'No match'
      : !ids.has(row.cluster) ? 'Unknown cluster' : counts.get(row.cluster) > 1 ? 'Cluster assigned more than once' : '';
  }
  const missing = clusters.filter(cluster => !counts.has(cluster.uid));
  return { rows, missing, ready: rows.every(row => !row.error) && missing.length === 0 };
}

export function resolveTree(tree, clusters, assignments = new Map()) {
  const matches = matchTreeLeaves(tree, clusters, assignments);
  if (!matches.ready) throw new Error('Match every tree leaf to a different cluster and include all plotted clusters.');
  const resolved = structuredClone(tree);
  const mapping = new Map(matches.rows.map(row => [row.uid, row.cluster]));
  for (const node of resolved.nodes) if (!node.children) node.cluster = mapping.get(node.uid);
  validateTree(resolved, new Set(clusters.map(cluster => cluster.uid)));
  return resolved;
}

function inputClusterKey(cluster, separateRecords) {
  const split = separateRecords && /\.(gb|gbk|genbank|gbf|gbff)$/i.test(cluster.source_file);
  return JSON.stringify([cluster.source_file, split ? cluster.loci.map(locus => locus.name) : null]);
}

// File/record identities survive changes in selection order and analysis IDs.
export function identifyInputClusters(clusters, separateRecords = false) {
  const identified = clusters.map(cluster => ({ ...cluster, uid: inputClusterKey(cluster, separateRecords) }));
  if (new Set(identified.map(cluster => cluster.uid)).size !== identified.length) {
    throw new Error('Tree matching requires distinct record identifiers within each split GenBank file.');
  }
  return identified;
}

export function inputTreeClusters(filesWithRecords, separateRecords = false) {
  const clusters = filesWithRecords.flatMap(({ file, records }) => {
    const base = { name: file.name.replace(/\.[^.]+$/, ''), source_file: file.name, loci: records };
    return separateRecords && /\.(gb|gbk|genbank|gbf|gbff)$/i.test(file.name)
      ? records.map(record => ({ ...base, name: record.organism || record.name, loci: [record] })) : [base];
  });
  return identifyInputClusters(clusters, separateRecords);
}

export function bindInputTree(tree, clusters, separateRecords = false) {
  const identified = identifyInputClusters(clusters, separateRecords);
  const ids = new Map(identified.map((cluster, i) => [cluster.uid, clusters[i].uid]));
  const bound = structuredClone(tree);
  for (const node of bound.nodes) {
    if (!node.children) {
      if (!ids.has(node.cluster)) throw new Error('Tree inputs changed during parsing. Resolve the tree matches again.');
      node.cluster = ids.get(node.cluster);
    }
  }
  validateTree(bound, new Set(clusters.map(cluster => cluster.uid)));
  return bound;
}
