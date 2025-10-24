import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('synthetic cell generation', () => {
  test('toggleCellVisibilityEnhanced synthesizes cell when missing', () => {
    const state = createMoleculeState({
      elements: ['C'],
      positions: [{ x: 2, y: 3, z: 4 }],
      cell: null,
    });
    expect(state.cell.enabled).toBe(false);
    state.toggleCellVisibilityEnhanced();
    expect(state.showCell).toBe(true);
    expect(state.cell.enabled).toBe(true);
    expect(state.cell.synthetic).toBe(true);
  });
});
