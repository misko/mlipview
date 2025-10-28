import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';

describe('selection highlight plumbing', () => {
  test('atom selection updates selection state', () => {
    const state = createMoleculeState({
      elements: ['O', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
    });
    const selection = createSelectionService(state);
    selection.clickAtom({ index: 0 });
    expect(selection.get()?.kind).toBe('atom');
    selection.clear();
    expect(selection.get()?.kind).toBeNull();
  });
});
