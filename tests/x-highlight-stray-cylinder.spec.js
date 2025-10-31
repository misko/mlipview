import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.ts';

describe('bond selection lifecycle', () => {
  test('selecting and clearing bond updates selection service', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.2, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const selection = createSelectionService(state);
    selection.clickBond({ i: 0, j: 1, key: '0-1', index: 0 });
    expect(selection.get()?.kind).toBe('bond');
    selection.clear();
    expect(selection.get()?.kind).toBeNull();
  });
});
