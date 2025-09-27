// public/forcefield/fairchem.js
// HuggingFace FAIR-Chem remote inference adapter.
// Lightweight JSON POST wrapper modeled after hf_calculator.py behavior.
// NOTE: Endpoint & auth token must be configured by user.

import { atomsToZ, atomsToXYZ, rawForcesToVectors, validateZXYZ } from './interface.js';

const DEFAULT_ENDPOINT = 'https://api-inference.huggingface.co/models/facebook/fairchem_mace_large';

function buildHeaders(token) {
  const h = { 'Content-Type': 'application/json' };
  if (token) h['Authorization'] = `Bearer ${token}`;
  return h;
}

async function postJSON(endpoint, payload, token) {
  const res = await fetch(endpoint, {
    method: 'POST',
    headers: buildHeaders(token),
    body: JSON.stringify(payload)
  });
  if (!res.ok) {
    let txt = '';
    try { txt = await res.text(); } catch {}
    throw new Error(`FAIR-Chem request failed: ${res.status} ${res.statusText} :: ${txt}`);
  }
  return res.json();
}

// Token handling: An optional file `public/hf-token.js` (ignored by git) can set window.HF_TOKEN = "hf_xxx".
// This avoids committing secrets. Example file content:
//   window.HF_TOKEN = 'hf_XXXXXXXXXXXXXXXXXXXXXXXXXXXX';
// Provide a token via that file or via an environment-driven injection mechanism during build/deploy.

export function createFairChemForceField({
  molecule = null,
  endpoint = DEFAULT_ENDPOINT,
  hfToken = (typeof window !== 'undefined' && window.HF_TOKEN)
    ? window.HF_TOKEN
    : null,
  responseMap = { energy: 'total_energy', forces: 'forces' },
  units = { energy: 'eV', length: 'Angstrom' },
  defaultCharge = 0,
  defaultSpin = 1
} = {}) {
  const boundAtoms = molecule ? molecule.atoms : null;

  async function computeRaw({ Z, xyz }) {
    validateZXYZ(Z, xyz);
    const payload = {
      atomic_numbers: Z,
      coordinates: xyz,
      charge: defaultCharge,
      spin_multiplicity: defaultSpin
    };
    const data = await postJSON(endpoint, payload, hfToken);
    const energy = data[responseMap.energy];
    const forces = data[responseMap.forces];
    if (!Array.isArray(forces) || forces.length !== Z.length) {
      throw new Error('Invalid forces array length from FAIR-Chem response');
    }
    return { energy, forces };
  }

  async function compute() {
    if (!boundAtoms) throw new Error('No bound molecule provided to FAIR-Chem force field');
    const Z = atomsToZ(boundAtoms);
    const xyz = atomsToXYZ(boundAtoms);
    const { energy, forces } = await computeRaw({ Z, xyz });
    return { energy, forces: rawForcesToVectors(forces) };
  }

  return {
    kind: 'fairchem',
    units,
    meta: { endpoint, provider: 'huggingface' },
    compute,
    computeRaw
  };
}
