import { elInfo } from "../elements.js";
import { buildAdjacency } from "../graph.js";

/**
 * Pure bond computation (no Babylon.js) returning bond list with opacity.
 * Inputs:
 *  atoms: [{ element, pos: BABYLON.Vector3 }]
 *  options:
 *    bondScale: cutoff multiplier for initial candidate bonds
 *    baseBondList: optional precomputed connectivity to seed ring detection [{i,j}]
 *    ringOptions: { aromaticScale, ringMin, ringMax }
 *    defaultFade: { min, max }   // fade region for non-ring bonds (in multiples of expected length)
 *    minOpacity: threshold below which bonds are discarded
 *    quantize: optional step to quantize opacity (e.g. 0.05); if 0 or falsy, no quantization
 *    includeMeta: if true include distance & expectedLength in each bond record
 * Returns:
 *  { bonds: [{ i,j,pairKey,opacity, ...(meta) }], ringFlags: boolean[] }
 */
export function computeBonds(atoms, options = {}) {
  const {
    bondScale = 1.15,
    baseBondList = null,
    ringOptions = { aromaticScale: 0.98, ringMin: 0.95, ringMax: 1.05 },
    defaultFade = { min: 0.9, max: 1.5 },
    minOpacity = 0.02,
    quantize = 0, // quantization step for opacity; 0 => disabled
    includeMeta = false
  } = options;

  const n = atoms.length;

  // Helper: element pair key (sorted)
  function pairKeyOf(i, j) {
    const a = atoms[i].element, b = atoms[j].element;
    return [a, b].sort().join("-");
  }

  // --- Phase 1: establish a base bond list for ring detection if not supplied.
  let seedBondList = baseBondList;
  if (!seedBondList) {
    seedBondList = [];
    for (let i = 0; i < n; i++) {
      const ai = atoms[i]; const infoi = elInfo(ai.element);
      for (let j = i + 1; j < n; j++) {
        const aj = atoms[j]; const infoj = elInfo(aj.element);
        const cutoff = (infoi.covRad + infoj.covRad) * bondScale;
        if (BABYLON.Vector3.Distance(ai.pos, aj.pos) <= cutoff) seedBondList.push({ i, j });
      }
    }
  }

  // --- Phase 2: ring detection (simple DFS cycles) using adjacency of seed bonds
  const adj = buildAdjacency(seedBondList, n);
  const ringFlags = findRingAtoms(adj);

  // --- Phase 3: compute opacity for every possible candidate pair (prune by large distance)
  const bonds = [];
  for (let i = 0; i < n; i++) {
    const ai = atoms[i]; const infoi = elInfo(ai.element);
    for (let j = i + 1; j < n; j++) {
      const aj = atoms[j]; const infoj = elInfo(aj.element);
      const expectedBase = (infoi.covRad + infoj.covRad) * bondScale;
      const distance = BABYLON.Vector3.Distance(ai.pos, aj.pos);

      // Quick rejection: outside an extended window (2.0x expected) no bond
      if (distance > expectedBase * 2.0) continue;

      const pairKey = pairKeyOf(i, j);
      let expectedBondLength = expectedBase;
      let minEdge, maxEdge; // absolute distance thresholds for zero opacity

      // Aromatic / ring C-C tighter window
      if (pairKey === 'C-C' && ringFlags[i] && ringFlags[j]) {
        const { aromaticScale, ringMin, ringMax } = ringOptions;
        expectedBondLength *= aromaticScale;
        minEdge = expectedBondLength * ringMin;
        maxEdge = expectedBondLength * ringMax;
      } else {
        minEdge = expectedBondLength * defaultFade.min;
        maxEdge = expectedBondLength * defaultFade.max;
      }

      let opacity = 1.0;
      if (distance < minEdge || distance > maxEdge) {
        opacity = 0.0;
      } else {
        if (distance < expectedBondLength) {
          const denom = expectedBondLength - minEdge;
          opacity = denom > 1e-12 ? (distance - minEdge) / denom : 1.0;
        } else if (distance > expectedBondLength) {
          const denom = maxEdge - expectedBondLength;
          opacity = denom > 1e-12 ? (maxEdge - distance) / denom : 1.0;
        } else {
          opacity = 1.0;
        }
      }

      if (opacity <= minOpacity) continue;
      if (quantize && quantize > 0) {
        const inv = 1 / quantize;
        opacity = Math.round(opacity * inv) / inv;
      }

      const rec = { i, j, pairKey, opacity };
      if (includeMeta) {
        rec.distance = distance;
        rec.expectedLength = expectedBondLength;
        rec.inRing = !!(ringFlags[i] && ringFlags[j]);
      }
      bonds.push(rec);
    }
  }

  return { bonds, ringFlags };
}

// Internal: simple cycle marking via DFS back-edge detection
function findRingAtoms(adj) {
  const n = adj.length;
  const inRing = new Array(n).fill(false);
  const visited = new Array(n).fill(false);

  function dfs(u, parent, stack, onStack) {
    visited[u] = true;
    stack.push(u);
    onStack[u] = stack.length - 1;
    for (const v of adj[u]) {
      if (v === parent) continue;
      if (!visited[v]) {
        dfs(v, u, stack, onStack);
      } else if (onStack[v] !== -1) {
        for (let k = onStack[v]; k < stack.length; k++) inRing[stack[k]] = true;
      }
    }
    onStack[u] = -1;
    stack.pop();
  }

  for (let i = 0; i < n; i++) {
    if (!visited[i]) {
      const stack = [];
      const onStack = new Array(n).fill(-1);
      dfs(i, -1, stack, onStack);
    }
  }
  return inRing;
}
