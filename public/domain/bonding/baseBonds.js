import { covalentRadius } from './covalent.js';
import {
  R_SCALE,
  SWITCH_IN,
  SWITCH_OUT,
  MAX_NEIGHBORS,
  RING_MIN,
  RING_MAX,
  RING_W_MIN,
  RING_BOOST,
  OPACITY_GAMMA,
  OPACITY_FLOOR,
} from './constants.js';
import { vsub, vlen, vcross, vadd, vscale, keyForPair } from './math.js';

const smoothstep = (a, b, x) => {
  const t = Math.min(1, Math.max(0, (x - a) / (b - a)));
  return t * t * (3 - 2 * t);
};

function bondWeight(distance, rc) {
  const r0 = SWITCH_IN * rc;
  const r1 = SWITCH_OUT * rc;
  return 1 - smoothstep(r0, r1, distance);
}

function planarity(indices, positions) {
  if (indices.length < 3) return 0;
  const pts = indices.map((idx) => positions[idx]);
  let centroid = [0, 0, 0];
  for (const p of pts) centroid = vadd(centroid, p);
  centroid = vscale(centroid, 1 / pts.length);
  let sumCross = [0, 0, 0];
  let denom = 0;
  for (let i = 0; i < pts.length; i++) {
    const a = vsub(pts[i], centroid);
    const b = vsub(pts[(i + 1) % pts.length], centroid);
    const cr = vcross(a, b);
    sumCross = vadd(sumCross, cr);
    denom += vlen(cr);
  }
  if (denom < 1e-9) return 0;
  return Math.min(1, vlen(sumCross) / denom);
}

function detectRings({ edges, nbrLookup, positions }) {
  const ringScoreByEdge = {};
  const seenCycleCanon = new Set();

  function dfs(start, current, path, visited, parent) {
    if (path.length > RING_MAX) return;
    const neighbors = nbrLookup.get(current) || [];
    for (const { nb, weight } of neighbors) {
      if (nb === parent) continue;
      const idxInPath = path.indexOf(nb);
      if (idxInPath !== -1) {
        const cycle = path.slice(idxInPath);
        if (cycle.length >= RING_MIN && cycle.length <= RING_MAX) {
          const canonical = cycle
            .slice()
            .sort((a, b) => a - b)
            .join(',');
          if (!seenCycleCanon.has(canonical)) {
            seenCycleCanon.add(canonical);
            let wsum = 0;
            for (let i = 0; i < cycle.length; i++) {
              const a = cycle[i];
              const b = cycle[(i + 1) % cycle.length];
              wsum += edges.get(keyForPair(a, b))?.weight ?? 0;
            }
            const wmean = wsum / cycle.length;
            const plan = planarity(cycle, positions);
            const score = wmean * plan;
            if (score > 0) {
              for (let i = 0; i < cycle.length; i++) {
                const a = cycle[i];
                const b = cycle[(i + 1) % cycle.length];
                const edgeKey = keyForPair(a, b);
                const prev = ringScoreByEdge[edgeKey] ?? 0;
                ringScoreByEdge[edgeKey] = Math.max(prev, score);
              }
            }
          }
        }
        continue;
      }
      if (visited.has(nb)) continue;
      if (nb < start) continue;
      visited.add(nb);
      dfs(start, nb, [...path, nb], visited, current);
      visited.delete(nb);
    }
  }

  const sortedKeys = [...nbrLookup.keys()].sort((a, b) => a - b);
  for (const carbonIdx of sortedKeys) {
    const visited = new Set([carbonIdx]);
    dfs(carbonIdx, carbonIdx, [carbonIdx], visited, -1);
  }
  return ringScoreByEdge;
}

export function computeBaseBonds(atoms, options = {}) {
  if (!Array.isArray(atoms) || !atoms.length) {
    return { bonds: [], diagnostics: { neighborCounts: [], prunedPairs: 0 } };
  }
  const debug = typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;
  const maxNeighbors = Number.isFinite(options.maxNeighbors)
    ? Math.max(1, options.maxNeighbors | 0)
    : MAX_NEIGHBORS;
  const positions = atoms.map((a) => a.pos);
  const elements = atoms.map((a) => a.element);
  const n = atoms.length;
  const neighbors = Array.from({ length: n }, () => []);
  let prunedPairs = 0;
  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const rc =
        R_SCALE * ((covalentRadius(elements[i]) ?? 0.8) + (covalentRadius(elements[j]) ?? 0.8));
      const distance = vlen(vsub(positions[i], positions[j]));
      if (distance < rc * SWITCH_OUT * 1.1) {
        neighbors[i].push(j);
        neighbors[j].push(i);
      }
    }
  }
  for (let i = 0; i < n; i++) {
    if (neighbors[i].length > maxNeighbors) {
      neighbors[i].sort((a, b) => {
        const da = vlen(vsub(positions[i], positions[a]));
        const db = vlen(vsub(positions[i], positions[b]));
        return da - db;
      });
      prunedPairs += neighbors[i].length - maxNeighbors;
      neighbors[i] = neighbors[i].slice(0, maxNeighbors);
      if (debug) {
        console.log('[Bonding][neighborTrim]', {
          atom: i,
          original: neighbors[i].length + (prunedPairs ? prunedPairs : 0),
          trimmedTo: maxNeighbors,
        });
      }
    }
  }

  const edgeMap = new Map();
  for (let i = 0; i < n; i++) {
    for (const j of neighbors[i]) {
      if (j <= i) continue;
      const rc =
        R_SCALE * ((covalentRadius(elements[i]) ?? 0.8) + (covalentRadius(elements[j]) ?? 0.8));
      const distance = vlen(vsub(positions[i], positions[j]));
      const weight = bondWeight(distance, rc);
      if (weight <= 0) continue;
      edgeMap.set(keyForPair(i, j), {
        i,
        j,
        distance,
        weight,
      });
    }
  }

  const ccAdjacency = new Map();
  for (let i = 0; i < n; i++) {
    if (elements[i] === 'C') {
      ccAdjacency.set(i, []);
    }
  }
  for (const edge of edgeMap.values()) {
    if (elements[edge.i] === 'C' && elements[edge.j] === 'C' && edge.weight >= RING_W_MIN) {
      ccAdjacency.get(edge.i)?.push({ nb: edge.j, weight: edge.weight });
      ccAdjacency.get(edge.j)?.push({ nb: edge.i, weight: edge.weight });
    }
  }

  const ringScores =
    ccAdjacency.size > 0
      ? detectRings({ edges: edgeMap, nbrLookup: ccAdjacency, positions })
      : {};

  const bonds = [];
  for (const edge of edgeMap.values()) {
    const key = keyForPair(edge.i, edge.j);
    const ringScore = ringScores[key] ?? 0;
    const boosted = edge.weight * (1 + RING_BOOST * ringScore);
    const opacity = Math.max(OPACITY_FLOOR, Math.min(1, Math.pow(boosted, OPACITY_GAMMA)));
    bonds.push({
      i: edge.i,
      j: edge.j,
      length: edge.distance,
      weight: edge.weight,
      opacity,
      inRing: ringScore > 0.05,
    });
  }

  if (debug) {
    const weights = bonds.map((b) => b.weight);
    console.log('[Bonding][baseComplete]', {
      atoms: n,
      bonds: bonds.length,
      prunedPairs,
      weightMin: weights.length ? Math.min(...weights) : null,
      weightMax: weights.length ? Math.max(...weights) : null,
    });
  }

  return {
    bonds,
    diagnostics: {
      neighborCounts: neighbors.map((list) => list.length),
      prunedPairs,
    },
  };
}
