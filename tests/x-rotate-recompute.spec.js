import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.ts';
import { createBondService } from '../public/domain/bondService.ts';
import { createSelectionService } from '../public/domain/selectionService.ts';

function buildLinearChain() {
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

describe('x-rotate-recompute', () => {
  test('rotateBond increments bond version counter', () => {
    const state = buildLinearChain();
    const bondService = createBondService(state);
    const manip = createManipulationService(state, { bondService });
    const selection = createSelectionService(state);

    const initialVersion = state.versions.bonds;
    selection.clickBond({ i: 0, j: 1, key: '0-1', index: 0 });
    manip.rotateBond(Math.PI / 6);
    expect(state.versions.bonds).toBeGreaterThan(initialVersion);
  });
});
