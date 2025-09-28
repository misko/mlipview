// public/state.js
// Simplified persistent molecule state (rotation history, rebuild/export removed).

import { buildAdjacency } from "./graph.js";

export function createStateStore(molecule) {
  const atomsAPI = molecule.atoms;   // [{type, index, pos, scale, mesh}]
  const n = atomsAPI.length;

  // Immutable element list
  const elements = atomsAPI.map(a => a.type);

  // Deep copy initial positions as plain vectors
  const initial = atomsAPI.map(a => a.pos.clone());

  // Bonds (indices) from the molecule builder
  const bonds = molecule.bonds.map(({ i, j }) => ({ i, j }));

  // Adjacency for side-set queries
  let adj = buildAdjacency(bonds, n);

  // Rotation history & advanced editing features removed (simplified state store).

  function getStateJSON() {
    return JSON.stringify({ elements, bonds }, null, 2);
  }

  function loadStateJSON(json) { /* no-op after simplification */ }

  function refreshAdjacency() {
    // Call this if you run molecule.recomputeBonds()
    adj = buildAdjacency(molecule.bonds, n);
  }

  return { molecule, elements, bonds, initial, getStateJSON, loadStateJSON, refreshAdjacency };
}
