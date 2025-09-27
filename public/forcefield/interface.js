// public/forcefield/interface.js
// Standard interface & helpers for force fields / MLIPs.
// Goal: unify classical (LJ) and remote ML (FAIR-Chem) backends.
//
// Contract (simplified):
//   Factory returns an object implementing:
//     kind: string                // identifier e.g. 'lj', 'fairchem'
//     units: { energy, length }   // optional metadata
//     meta:  { ... }              // backend specific metadata
//     supportsBatch?: boolean
//     compute(): { energy:number, forces: BABYLON.Vector3[] }
//       - Uses bound molecule supplied at factory time
//     computeRaw({ Z:int[], xyz:number[][] }): Promise<{ energy:number, forces:number[][] }>
//       - Pure numeric form; no Babylon dependencies (except via caller converting)
//
// Notes:
//  - atomic numbers Z follow standard periodic table mapping below.
//  - xyz units are Angstrom unless noted otherwise.

export const ELEMENT_Z = {
  H: 1, He: 2,
  C: 6, N: 7, O: 8, F: 9, Ne: 10,
  S: 16,
};

export const Z_ELEMENT = Object.fromEntries(Object.entries(ELEMENT_Z).map(([k,v]) => [v,k]));

export function symbolToZ(sym) {
  return ELEMENT_Z[sym] || 0;
}

export function zToSymbol(z) {
  return Z_ELEMENT[z] || 'X';
}

export function atomsToZ(atoms) {
  return atoms.map(a => symbolToZ(a.type || a.element || a.symbol));
}

export function atomsToXYZ(atoms) {
  return atoms.map(a => [a.pos.x, a.pos.y, a.pos.z]);
}

export function rawForcesToVectors(forces) {
  return forces.map(f => new BABYLON.Vector3(f[0], f[1], f[2]));
}

export function validateZXYZ(Z, xyz) {
  if (!Array.isArray(Z) || !Array.isArray(xyz)) throw new Error('Z and xyz must be arrays');
  if (Z.length !== xyz.length) throw new Error('Z and xyz length mismatch');
  for (let i=0;i<xyz.length;i++) {
    const r = xyz[i];
    if (!Array.isArray(r) || r.length !== 3 || !r.every(n => Number.isFinite(n))) {
      throw new Error('Invalid coordinate at index '+i);
    }
  }
}
