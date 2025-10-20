import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

function simpleState() {
  // Linear triatom chain 0-1-2 so rotating around bond 0-1 affects atom 0 or 1 side
  return createMoleculeState({
    elements: ['C', 'C', 'H'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
      { x: 2, y: 0, z: 0 },
    ],
    bonds: [
      { i: 0, j: 1 },
      { i: 1, j: 2 },
    ],
  });
}

describe('rotateBond recompute', () => {
  test('recomputes bonds after rotation', () => {
    const st = simpleState();
    const bondService = createBondService(st);
    const manip = createManipulationService(st, { bondService });
    // Simulate selecting bond 0-1 with orientation metadata
    st.selection = { kind: 'bond', data: { i: 0, j: 1, orientation: 'i->j' } };
    const before = st.versions.bonds;
    manip.rotateBond(Math.PI / 4);
    const after = st.versions.bonds;
    expect(after).toBeGreaterThan(before);
  });
});
