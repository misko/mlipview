import { computeBaseBonds } from './baseBonds.ts';
import { augmentWithPeriodicImages } from './periodicAugment.ts';
import type {
  BaseBondResult,
  BondRecord,
  BondComputationOptions,
  ComputeBondsResult,
  NormalizedAtom,
  PeriodicAugmentResult,
  PeriodicCell,
  Vec3,
} from './types.ts';

function normalizeAtoms(elements: string[], positions: Array<unknown>): NormalizedAtom[] {
  const atoms: NormalizedAtom[] = [];
  const len = Math.min(elements.length, positions.length);
  for (let i = 0; i < len; i++) {
    const el = elements[i];
    const pos = positions[i];
    if (!Array.isArray(pos) || pos.length < 3) continue;
    const x = Number(pos[0]) || 0;
    const y = Number(pos[1]) || 0;
    const z = Number(pos[2]) || 0;
    atoms.push({ element: el, pos: [x, y, z] });
  }
  return atoms;
}

interface ComputeBondsInput {
  elements?: string[];
  positions?: Array<unknown>;
  cell?: PeriodicCell | null;
  options?: BondComputationOptions;
}

export function computeBonds({
  elements = [],
  positions = [],
  cell = null,
  options = {},
}: ComputeBondsInput = {}): ComputeBondsResult {
  const atoms = normalizeAtoms(elements, positions);
  if (!atoms.length) {
    return {
      bonds: [],
      ghostAtoms: [],
      ghostBondMeta: [],
      diagnostics: {},
    };
  }

  const base = computeBaseBonds(atoms, options) as BaseBondResult;
  const baseBonds: BondRecord[] = base.bonds.map((bond) => ({
    i: bond.i,
    j: bond.j,
    length: bond.length ?? null,
    weight: bond.weight ?? null,
    opacity: bond.opacity,
    inRing: bond.inRing,
    crossing: false,
    imageDelta: [0, 0, 0] as Vec3,
    cellOffsetA: [0, 0, 0] as Vec3,
    cellOffsetB: [0, 0, 0] as Vec3,
  }));

  if (!cell || !cell.enabled) {
    return {
      bonds: baseBonds,
      ghostAtoms: [],
      ghostBondMeta: [],
      diagnostics: { base: base.diagnostics },
    };
  }

  const periodic = augmentWithPeriodicImages({
    atoms,
    cell,
    options,
  }) as PeriodicAugmentResult;

  const bonds = periodic.bonds.length ? periodic.bonds : baseBonds;

  if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true) {
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

declare global {
  interface Window {
    __MLIPVIEW_DEBUG_BONDS?: boolean;
  }
}
