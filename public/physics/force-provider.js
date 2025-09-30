// Unified force provider abstraction (local LJ + FairChem remote)
// Exposes a consistent async compute interface returning energy, forces, optional stress.
// Units: energy eV, forces eV/Å, stress eV/Å^3.

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
    // Ensure force array exists (mirrors integrators initialization logic)
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

// FairChem remote provider hitting FastAPI endpoint. By default we assume the backend is exposed
// on the SAME origin (served or reverse-proxied) and we call a relative path '/simple_calculate'.
// You can still override with a full baseUrl (e.g. http://host:8000) by passing { baseUrl } or
// setting window.__FAIRCHEM_URL. If baseUrl is omitted we use relative fetch for zero CORS hassle.
export function createFairChemProvider({ baseUrl } = {}) {
  if (!baseUrl && typeof window !== 'undefined' && window.__FAIRCHEM_URL) {
    baseUrl = window.__FAIRCHEM_URL; // explicit global override
  }
  const capabilities = { supportsStress: true, supportsCellRelax: true, name: 'fairchem-remote' };
  async function compute({ elements = [], positions = [], cell = null, wantStress = true } = {}) {
    const Z = ensureZ(elements);
    const coords = positions.map(p => Array.isArray(p) ? p : [p.x, p.y, p.z]);
    const props = ['energy', 'forces'];
    if (wantStress) props.push('stress');
    const body = { atomic_numbers: Z, coordinates: coords, properties: props };
    if (cell) body.cell = cell; // Expect 3x3 array
    let respData = null;
    try {
      const url = baseUrl ? (baseUrl.replace(/\/$/, '') + '/simple_calculate') : '/simple_calculate';
      const resp = await fetch(url, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(body)
      });
      if (!resp.ok) throw new Error('HTTP ' + resp.status);
      respData = await resp.json();
    } catch (e) {
      // Network or server failure -> provide deterministic mock
      const n = Z.length;
      return { energy: n * 0.05, forces: Z.map(()=>[0,0,0]), stress: wantStress ? { tensor:null, voigt:null } : null, error: e.message };
    }
    const results = respData.results || respData || {};
    const energy = results.energy;
    const forces = (results.forces || []).map(f => [f[0], f[1], f[2]]);
    let stress = null;
    if (wantStress && Array.isArray(results.stress)) {
      const v = results.stress;
      if (v.length === 6) {
        stress = { tensor: { xx:v[0], yy:v[1], zz:v[2], yz:v[3], xz:v[4], xy:v[5] }, voigt: v.slice() };
      }
    }
    return { energy, forces, stress };
  }
  return { compute, capabilities };
}

// Generic factory based on type string.
export function createForceProvider(kind, opts) {
  if (kind === 'local') return createLocalLJProvider(opts.molState, opts);
  if (kind === 'fairchem') return createFairChemProvider(opts);
  throw new Error('Unknown force provider kind: ' + kind);
}
