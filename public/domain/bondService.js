import { computeBonds } from './bonding/index.js';
import { __count } from '../util/funcCount.js';

export function createBondService(molState) {
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
    // Store opacity for renderer so it can apply group alpha (harmless to existing logic using only i,j)
    // Preserve all bonds with per-instance opacity; renderer will hide crossing ones visually.
    const DBG =
      typeof window !== 'undefined' &&
      (window.__MLIP_DEBUG_STRETCH === true || /[?&]bondStretchDebug=1/.test(window.location?.search || ''));
    if (DBG) {
      try {
        let min = Number.POSITIVE_INFINITY;
        let max = Number.NEGATIVE_INFINITY;
        let softish = 0;
        const samples = [];
        for (const b of bonds) {
          const op = typeof b.opacity === 'number' ? b.opacity : 1;
          if (op < min) min = op;
          if (op > max) max = op;
          if (op < 0.99) softish++;
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
    const DEBUG_BONDS =
      typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_BONDS === true;
    if (DEBUG_BONDS) {
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
    const ghostDebug =
      typeof window !== 'undefined' && window.__MLIP_DEBUG_GHOST_BONDS === true;

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
      ? result.ghostAtoms.map((ghost) => ({
          atomIndex: ghost.atomIndex,
          shift: ghost.shift ? ghost.shift.slice() : [0, 0, 0],
          position: Array.isArray(ghost.position)
            ? ghost.position.slice(0, 3)
            : [ghost.position?.x || 0, ghost.position?.y || 0, ghost.position?.z || 0],
        }))
      : [];
    molState.ghostBondMeta = Array.isArray(result?.ghostBondMeta)
      ? result.ghostBondMeta.map((meta) => ({
          base: { i: meta.base?.i ?? 0, j: meta.base?.j ?? 0 },
          shiftA: Array.isArray(meta.shiftA) ? meta.shiftA.slice(0, 3) : [0, 0, 0],
          shiftB: Array.isArray(meta.shiftB) ? meta.shiftB.slice(0, 3) : [0, 0, 0],
          imageDelta: Array.isArray(meta.imageDelta) ? meta.imageDelta.slice(0, 3) : [0, 0, 0],
        }))
      : [];

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
    return bonds;
  }
  function computePeriodicBonds() {
    __count('bondService#computePeriodicBonds');
    return recomputeAndStore();
  }
  return { computePeriodicBonds, recomputeAndStore };
}
