// Simplified local force provider (deprecated multi-provider abstraction removed).
import { createForceField } from './forcefield.js';

// Minimal periodic table subset for Z lookup (extend as needed)
const ELEMENT_Z = { H:1, C:6, N:7, O:8, S:16 };

function ensureZ(elements) {
  return elements.map(e => {
    if (typeof e === 'number') return e;
    const z = ELEMENT_Z[e];
    if (!z) throw new Error('Unknown element: ' + e);
    return z;
  });
}

export function createLocalLJProvider(molState, opts={}) {
  // We reuse existing local forcefield implementation.
  const ff = createForceField(molState, opts);
  const capabilities = { supportsStress: true, supportsCellRelax: false, name: 'local-lj' };
  async function compute({ wantStress = true, positions } = {}) {
    // If caller supplies positions array, sync into molState first
    if (positions && Array.isArray(positions) && positions.length === molState.positions.length) {
      for (let i=0;i<positions.length;i++) {
        const p = positions[i];
        const mp = molState.positions[i];
        mp.x = p[0]; mp.y = p[1]; mp.z = p[2];
      }
    }
  // Ensure force array exists
    if (!molState.dynamics.forces || !molState.dynamics.forces.length) {
      molState.dynamics.forces = molState.positions.map(()=>({x:0,y:0,z:0}));
    }
    ff.computeForces();
    // Extract forces as flat JS arrays to decouple from mutable objects.
    const forces = molState.dynamics.forces.map(f => [f.x, f.y, f.z]);
    let stress = null;
    if (wantStress && molState.dynamics.stress) {
      const s = molState.dynamics.stress;
      stress = { tensor: { xx:s.xx, yy:s.yy, zz:s.zz, xy:s.xy, xz:s.xz, yz:s.yz }, voigt: [s.xx, s.yy, s.zz, s.yz, s.xz, s.xy] };
    }
    return { energy: undefined, // local model presently accumulates energy but not exposed externally yet
      forces, stress };
  }
  return { compute, capabilities };
}

// NOTE: createFairChemProvider & createForceProvider removed as unused after cleanup.
