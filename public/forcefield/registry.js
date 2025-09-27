// public/forcefield/registry.js
// Simple registry to construct force fields by name.

import { createLJForceField } from '../physics_mock.js';
import { createFairChemForceField } from './fairchem.js';

const FACTORIES = {
  lj: ({ molecule }) => createLJForceField(molecule),
  fairchem: (opts) => createFairChemForceField(opts)
};

export function listForceFields() {
  return Object.keys(FACTORIES);
}

export function createForceField(kind, opts) {
  const f = FACTORIES[kind];
  if (!f) throw new Error(`Unknown force field kind: ${kind}`);
  return f(opts || {});
}
