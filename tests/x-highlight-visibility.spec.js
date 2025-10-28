import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('selection version tracking', () => {
  test('markSelectionChanged bumps version counter', () => {
    const state = createMoleculeState({
      elements: ['H'],
      positions: [{ x: 0, y: 0, z: 0 }],
    });
    const v0 = state.versions.selection;
    state.markSelectionChanged();
    expect(state.versions.selection).toBe(v0 + 1);
  });
});
