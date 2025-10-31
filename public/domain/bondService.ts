import { computeBonds } from './bonding/index.js';
import { __count } from '../util/funcCount.ts';
import type { ComputeBondsResult, GhostAtom, GhostBondMeta, PeriodicCell } from './bonding/types.ts';

interface PositionVec {
  x: number;
  y: number;
  z: number;
}

interface MoleculeState {
  elements: string[];
  positions: PositionVec[];
  cell?: PeriodicCell | null;
  bonds: Array<{
    i: number;
    j: number;
    length: number | null | undefined;
    weight: number | null | undefined;
    opacity: number;
    inRing?: boolean;
    crossing?: boolean;
    imageDelta?: number[];
    cellOffsetA?: number[];
    cellOffsetB?: number[];
  }>;
  ghostImages: Array<{
    atomIndex: number;
    shift: number[];
    position: number[];
  }>;
  ghostBondMeta: Array<{
    base: { i: number; j: number };
    shiftA: number[];
    shiftB: number[];
    imageDelta: number[];
  }>;
  showCell?: boolean;
  showGhostCells?: boolean;
  markBondsChanged: () => void;
}

function normalizeGhostAtom(atom: GhostAtom) {
  const shift = Array.isArray(atom.shift) ? atom.shift.slice(0, 3) : [0, 0, 0];
  const position = Array.isArray(atom.position)
    ? atom.position.slice(0, 3)
    : [
        (atom.position as { x?: number })?.x || 0,
        (atom.position as { y?: number })?.y || 0,
        (atom.position as { z?: number })?.z || 0,
      ];
  return { atomIndex: atom.atomIndex, shift, position };
}

function normalizeGhostMeta(meta: GhostBondMeta) {
  return {
    base: { i: meta.base?.i ?? 0, j: meta.base?.j ?? 0 },
    shiftA: Array.isArray(meta.shiftA) ? meta.shiftA.slice(0, 3) : [0, 0, 0],
    shiftB: Array.isArray(meta.shiftB) ? meta.shiftB.slice(0, 3) : [0, 0, 0],
    imageDelta: Array.isArray(meta.imageDelta) ? meta.imageDelta.slice(0, 3) : [0, 0, 0],
  };
}

function normalizeBondResult(result: ComputeBondsResult, molState: MoleculeState) {
  const bonds = Array.isArray(result?.bonds) ? result.bonds : [];
  molState.bonds = bonds.map((b) => ({
    i: b.i,
    j: b.j,
    length: b.length,
    weight: typeof b.weight === 'number' ? b.weight : null,
    opacity: typeof b.opacity === 'number' ? b.opacity : 1,
    inRing: !!b.inRing,
    crossing: !!b.crossing,
    imageDelta: Array.isArray(b.imageDelta) ? b.imageDelta.slice(0, 3) : [0, 0, 0],
    cellOffsetA: Array.isArray(b.cellOffsetA) ? b.cellOffsetA.slice(0, 3) : [0, 0, 0],
    cellOffsetB: Array.isArray(b.cellOffsetB) ? b.cellOffsetB.slice(0, 3) : [0, 0, 0],
  }));

  molState.ghostImages = Array.isArray(result?.ghostAtoms)
    ? result.ghostAtoms.map(normalizeGhostAtom)
    : [];

  molState.ghostBondMeta = Array.isArray(result?.ghostBondMeta)
    ? result.ghostBondMeta.map(normalizeGhostMeta)
    : [];
}

export function createBondService(molState: MoleculeState) {
  __count('bondService#createBondService');

  function recomputeAndStore() {
    __count('bondService#recomputeAndStore');
    const result = computeBonds({
      elements: molState.elements,
      positions: molState.positions.map((p) => [p.x, p.y, p.z]),
      cell: molState.cell,
      options: {},
    });

    const bonds = Array.isArray(result?.bonds) ? result.bonds : [];

    const stretchDebug =
      typeof window !== 'undefined' &&
      (window.__MLIP_DEBUG_STRETCH === true ||
        /[?&]bondStretchDebug=1/.test(window.location?.search || ''));

    if (stretchDebug) {
      try {
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        let softish = 0;
        const samples: Array<{ i: number; j: number; opacity: number | null; crossing: boolean }> = [];
        for (const b of bonds) {
          const op = typeof b.opacity === 'number' ? b.opacity : 1;
          if (op < min) min = op;
          if (op > max) max = op;
          if (op < 0.99) softish += 1;
        }
        for (let i = 0; i < Math.min(6, bonds.length); i++) {
          const b = bonds[i];
          samples.push({
            i: b.i,
            j: b.j,
            opacity: typeof b.opacity === 'number' ? Number(b.opacity.toFixed(4)) : null,
            crossing: !!b.crossing,
          });
        }
        console.log('[BondService][recompute]', {
          count: bonds.length,
          minOpacity: Number.isFinite(min) ? Number(min.toFixed(4)) : null,
          maxOpacity: Number.isFinite(max) ? Number(max.toFixed(4)) : null,
          belowThreshold: softish,
          sample: samples,
        });
      } catch (err) {
        console.warn('[BondService][recompute] debug log error', err);
      }
    }

    const bondsDebug = typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;
    if (bondsDebug) {
      const periodicDiag = result?.diagnostics?.periodic || {};
      console.log('[BondService][diagnostics]', {
        bonds: bonds.length,
        ghostAtoms: Array.isArray(result?.ghostAtoms) ? result.ghostAtoms.length : 0,
        ghostBondMeta: Array.isArray(result?.ghostBondMeta) ? result.ghostBondMeta.length : 0,
        base: result?.diagnostics?.base,
        periodic: periodicDiag,
      });
      if (Array.isArray(result?.ghostBondMeta) && result.ghostBondMeta.length) {
        console.log('[BondService][ghostMetaSample]', result.ghostBondMeta.slice(0, 8));
      }
    }

    normalizeBondResult(result, molState);

    const ghostDebug = typeof window !== 'undefined' && window.__MLIP_DEBUG_GHOST_BONDS === true;
    if (ghostDebug) {
      try {
        console.log('[BondService][ghostSummary]', {
          reason: 'recomputeAndStore',
          bonds: molState.bonds.length,
          ghostAtoms: molState.ghostImages.length,
          ghostBonds: molState.ghostBondMeta.length,
          showCell: !!molState.showCell,
          showGhostCells: !!molState.showGhostCells,
        });
      } catch (err) {
        console.warn('[BondService][ghostSummary] log failed', err);
      }
    }

    molState.markBondsChanged();
    return molState.bonds;
  }

  function computePeriodicBonds() {
    __count('bondService#computePeriodicBonds');
    return recomputeAndStore();
  }

  return { computePeriodicBonds, recomputeAndStore };
}

declare global {
  interface Window {
    __MLIP_DEBUG_STRETCH?: boolean;
    __MLIPVIEW_DEBUG_BONDS?: boolean;
    __MLIP_DEBUG_GHOST_BONDS?: boolean;
  }
}
