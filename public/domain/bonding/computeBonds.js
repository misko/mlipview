import { computeBaseBonds } from './baseBonds.js';
import { augmentWithPeriodicImages } from './periodicAugment.js';

function normalizeAtoms(elements, positions) {
  const atoms = [];
  for (let i = 0; i < elements.length; i++) {
    const el = elements[i];
    const pos = positions[i];
    if (!pos || !Array.isArray(pos) || pos.length < 3) continue;
    atoms.push({
      element: el,
      pos: [Number(pos[0]) || 0, Number(pos[1]) || 0, Number(pos[2]) || 0],
    });
  }
  return atoms;
}

export function computeBonds({ elements = [], positions = [], cell = null, options = {} } = {}) {
  const atoms = normalizeAtoms(elements, positions);
  if (!atoms.length) {
    return { bonds: [], ghostAtoms: [], ghostBondMeta: [], diagnostics: {} };
  }
  const base = computeBaseBonds(atoms, options);
  const baseBonds = base.bonds.map((bond) => ({
    ...bond,
    crossing: false,
    imageDelta: [0, 0, 0],
    cellOffsetA: [0, 0, 0],
    cellOffsetB: [0, 0, 0],
  }));

  if (!cell || !cell.enabled) {
    return {
      bonds: baseBonds,
      ghostAtoms: [],
      ghostBondMeta: [],
      diagnostics: { base: base.diagnostics },
    };
  }

  const periodic = augmentWithPeriodicImages({ atoms, cell, options });
  const bonds = periodic.bonds.length ? periodic.bonds : baseBonds;

  const debug = typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;
  if (debug) {
    console.log('[Bonding][computeSummary]', {
      baseCount: baseBonds.length,
      periodicCount: periodic.bonds.length,
      ghostAtoms: periodic.ghostAtoms.length,
      ghostBondMeta: periodic.ghostBondMeta.length,
    });
    if (!periodic.ghostBondMeta.length) {
      console.log('[Bonding][computeSummary] no ghostBondMeta generated', {
        elements: elements.length,
        cell,
      });
    }
  }

  return {
    bonds,
    ghostAtoms: periodic.ghostAtoms,
    ghostBondMeta: periodic.ghostBondMeta,
    diagnostics: {
      base: base.diagnostics,
      periodic: periodic.diagnostics,
    },
  };
}
