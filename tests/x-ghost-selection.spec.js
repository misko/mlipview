import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';

describe('selection metadata', () => {
  test('bond selection preserves orientation metadata', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 2, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const selection = createSelectionService(state);
    selection.clickBond({ i: 0, j: 1, key: '0-1', index: 0, orientation: 'ij' });
    const sel = selection.get();
    expect(sel?.kind).toBe('bond');
    expect(sel?.data.i).toBe(0);
    expect(sel?.data.j).toBe(1);
  });
});
