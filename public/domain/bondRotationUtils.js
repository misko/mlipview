import { orientationToSide } from '../selection-model.js';

export const ROTATION_BOND_OPACITY_MIN = 0.85;

export function computeBondRotationGroup(state, { i, j, orientation }) {
  if (!state || !Array.isArray(state.positions) || !Array.isArray(state.bonds)) {
    throw new Error('invalid state for bond rotation');
  }
  const side = orientationToSide(orientation);
  const anchor = side === 'i' ? j : i;
  const movingRoot = side === 'i' ? i : j;

  const adj = Array.from({ length: state.positions.length }, () => []);
  for (const b of state.bonds) {
    if (typeof b?.i !== 'number' || typeof b?.j !== 'number') continue;
    const opacity = b.opacity == null ? 1 : b.opacity;
    const isSelected = (b.i === i && b.j === j) || (b.i === j && b.j === i);
    if (isSelected || opacity >= ROTATION_BOND_OPACITY_MIN) {
      adj[b.i].push(b.j);
      adj[b.j].push(b.i);
    }
  }

  const visited = new Set([anchor]);
  const queue = [movingRoot];
  visited.add(movingRoot);
  const sideAtoms = [];
  while (queue.length) {
    const atom = queue.shift();
    sideAtoms.push(atom);
    const neighbors = adj[atom] || [];
    for (const nb of neighbors) {
      if (!visited.has(nb)) {
        visited.add(nb);
        queue.push(nb);
      }
    }
  }

  return {
    i,
    j,
    side,
    orientation,
    anchor,
    movingRoot,
    sideAtoms,
  };
}
