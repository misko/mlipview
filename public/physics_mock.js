// public/physics_mock.js
// Mock MLIP: LJ nonbonded + harmonic bonds.
// Usage: const mlip = createMockMLIP(molecule); const { energy, forces } = mlip.compute();

function pairKey(a, b) {
  return [a, b].sort().join("-");
}

// Some gentle, visualization-friendly parameters (units arbitrary but consistent)
const LJ_PARAMS = {
  "H-H": { eps: 0.02, sig: 2.2 },
  "C-H": { eps: 0.05, sig: 3.0 },
  "C-C": { eps: 0.08, sig: 3.4 },
  "C-N": { eps: 0.07, sig: 3.3 },
  "C-O": { eps: 0.09, sig: 3.1 },
  "C-S": { eps: 0.10, sig: 3.6 },
  "H-N": { eps: 0.04, sig: 2.9 },
  "H-O": { eps: 0.05, sig: 2.7 },
  "N-O": { eps: 0.07, sig: 3.0 },
  "N-N": { eps: 0.06, sig: 3.2 },
  "O-O": { eps: 0.10, sig: 3.0 },
  "S-H": { eps: 0.05, sig: 3.2 },
  "S-N": { eps: 0.08, sig: 3.6 },
  "S-O": { eps: 0.09, sig: 3.5 },
};

// Fallback if pair missing
const LJ_DEFAULT = { eps: 0.06, sig: 3.2 };

// Soft spring on bonds (kept low so it doesn't “fight” you while dragging)
const BOND_K = 0.15; // spring constant
const BOND_MAX = 6.0; // ignore silly long bonds

export function createMockMLIP(molecule) {
  const atoms = molecule.atoms; // [{type, pos, ...}] order = global index
  const n = atoms.length;
  
  // Performance tracking
  let computeCount = 0;

  // Map initial bond length r0 for harmonic term
  const r0 = new Map(); // key "i-j" (sorted) -> number
  for (const { i, j } of molecule.bonds) {
    const a = atoms[i].pos, b = atoms[j].pos;
    const key = i < j ? `${i}-${j}` : `${j}-${i}`;
    r0.set(key, BABYLON.Vector3.Distance(a, b));
  }

  function ljParams(a, b) {
    const k = pairKey(a, b);
    return LJ_PARAMS[k] || LJ_DEFAULT;
  }

  function compute() {
    computeCount++;
    console.log(`[MLIP DEBUG] compute() called (#${computeCount}) - this should only happen when atoms move or bonds rotate!`);
    
    const forces = Array.from({ length: n }, () => new BABYLON.Vector3(0, 0, 0));
    let E = 0.0;

    // --- Non-bonded LJ (naive O(N^2), fine for small molecules)
    for (let i = 0; i < n; i++) {
      const ai = atoms[i];
      for (let j = i + 1; j < n; j++) {
        const aj = atoms[j];
        const rij = aj.pos.subtract(ai.pos);
        const r2 = rij.lengthSquared();
        if (r2 < 1e-10) continue; // skip zero
        const r = Math.sqrt(r2);

        const { eps, sig } = ljParams(ai.type, aj.type);
        const sr = sig / r;
        const sr2 = sr * sr;
        const sr6 = sr2 * sr2 * sr2;
        const sr12 = sr6 * sr6;

        // Energy: 4eps (sr12 - sr6)
        const e = 4 * eps * (sr12 - sr6);
        E += e;

        // Force magnitude: d/dr of LJ -> 24eps/r (2 sr12 - sr6)
        const fmag = (24 * eps / r) * (2 * sr12 - sr6);
        const fij = rij.scale(fmag / r); // unit * fmag

        forces[i].subtractInPlace(fij);
        forces[j].addInPlace(fij);
      }
    }

    // --- Bond harmonic springs
    for (const { i, j } of molecule.bonds) {
      const key = i < j ? `${i}-${j}` : `${j}-${i}`;
      const rest = r0.get(key);
      const ri = atoms[i].pos, rj = atoms[j].pos;
      const d = rj.subtract(ri);
      const L = d.length();
      if (L < 1e-8 || L > BOND_MAX) continue;

      const x = L - rest;
      const eSpring = 0.5 * BOND_K * x * x;
      E += eSpring;

      const fmag = (BOND_K * x) / L;
      const f = d.scale(fmag); // direction from i->j

      forces[i].addInPlace(f);
      forces[j].subtractInPlace(f);
    }

    return { energy: E, forces };
  }

  return { compute };
}

// ---- New unified interface adapter (LJ + harmonic) ----
// Provides compute() (Babylon vectors) and computeRaw({Z,xyz}) (plain numbers)
// to match the generic MLIP/forcefield contract.
// Note: computeRaw currently omits harmonic bond term for simplicity; can be added if needed.
export function createLJForceField(molecule) {
  const base = createMockMLIP(molecule); // reuse existing implementation for scene-bound compute
  const atoms = molecule.atoms;
  const n = atoms.length;

  // Minimal element->Z map (extend as necessary)
  const ZMAP = { H:1, C:6, N:7, O:8, S:16 };
  const REV = Object.fromEntries(Object.entries(ZMAP).map(([k,v])=>[v,k]));

  function ljParams(a, b) {
    const k = pairKey(a, b);
    return LJ_PARAMS[k] || LJ_DEFAULT;
  }

  function computeRaw({ Z, xyz }) {
    if (!Z || !xyz) {
      // derive from bound atoms if not provided
      Z = atoms.map(a => ZMAP[a.type] || 0);
      xyz = atoms.map(a => [a.pos.x, a.pos.y, a.pos.z]);
    }
    if (Z.length !== xyz.length) throw new Error('Z/xyz length mismatch');
    const types = Z.map(z => REV[z] || 'X');
    const forces = Array.from({ length: Z.length }, () => [0,0,0]);
    let E = 0;
    for (let i=0;i<Z.length;i++) {
      const ti = types[i];
      const xi = xyz[i];
      for (let j=i+1;j<Z.length;j++) {
        const tj = types[j];
        const xj = xyz[j];
        const dx = xj[0]-xi[0];
        const dy = xj[1]-xi[1];
        const dz = xj[2]-xi[2];
        const r2 = dx*dx+dy*dy+dz*dz;
        if (r2 < 1e-12) continue;
        const r = Math.sqrt(r2);
        const { eps, sig } = ljParams(ti, tj);
        const sr = sig / r;
        const sr2 = sr*sr;
        const sr6 = sr2*sr2*sr2;
        const sr12 = sr6*sr6;
        const e = 4*eps*(sr12 - sr6);
        E += e;
        const fmag = (24*eps / r) * (2*sr12 - sr6);
        const fx = fmag * dx / r;
        const fy = fmag * dy / r;
        const fz = fmag * dz / r;
        forces[i][0] -= fx; forces[i][1] -= fy; forces[i][2] -= fz;
        forces[j][0] += fx; forces[j][1] += fy; forces[j][2] += fz;
      }
    }
    return { energy: E, forces };
  }

  function compute() { return base.compute(); }

  return {
    kind: 'lj',
    units: { energy: 'arb', length: 'arb' },
    meta: { source: 'local-lj' },
    compute,
    computeRaw
  };
}
