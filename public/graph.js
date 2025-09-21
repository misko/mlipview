// Build and use an undirected bond graph (adjacency list)

export function buildAdjacency(bonds, atomCount) {
  const adj = Array.from({ length: atomCount }, () => []);
  for (const { i, j } of bonds) {
    if (i < 0 || j < 0 || i >= atomCount || j >= atomCount) continue;
    adj[i].push(j);
    adj[j].push(i);
  }
  return adj;
}

/**
 * Return the set of atoms on the "side" of j when rotating bond (i–j),
 * i.e., all atoms reachable from j WITHOUT crossing i.
 */
export function sideAtoms(adj, i, j) {
  const visited = new Array(adj.length).fill(false);
  const side = new Set();
  const q = [j];
  visited[i] = true; // forbid passing through i
  visited[j] = true;
  side.add(j);

  while (q.length) {
    const u = q.shift();
    for (const v of adj[u]) {
      if (!visited[v]) {
        visited[v] = true;
        side.add(v);
        q.push(v);
      }
    }
  }
  return side;
}
