/**
 * This test is a light-weight simulation that ensures calling manipulation.rotateBond updates positions version.
 * Full DOM interaction (clicking + / - buttons) would require jsdom layout; here we invoke manipulation directly
 * after selecting a bond to assert state mutation.
 */
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

describe('bond rotation manipulation', () => {
  test('rotateBond increments positions version', () => {
    const state = createMoleculeState({
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
    const sel = createSelectionService(state);
    const manip = createManipulationService(state);
    sel.clickBond({ i: 0, j: 1, key: 'C-C', index: 0 });
    const before = state.versions.positions;
    const changed = manip.rotateBond(0.2);
    expect(changed).toBe(true);
    expect(state.versions.positions).toBe(before + 1);
  });
});
