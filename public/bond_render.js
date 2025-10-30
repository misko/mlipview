import { computeBaseBonds } from './domain/bonding/baseBonds.js';

export function computeBondsNoState(atoms) {
  if (!Array.isArray(atoms) || !atoms.length) return [];
  const normalized = atoms
    .map((atom) => ({
      element: atom.element,
      pos: Array.isArray(atom.pos)
        ? atom.pos.slice(0, 3).map((v) => Number(v) || 0)
        : [0, 0, 0],
    }))
    .filter((atom) => !!atom.element);
  const result = computeBaseBonds(normalized);
  return Array.isArray(result?.bonds) ? result.bonds : [];
}
