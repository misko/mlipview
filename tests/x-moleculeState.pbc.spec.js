import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('x-moleculeState-pbc', () => {
  test('toggleCellVisibilityEnhanced synthesizes padded cell and wraps atoms', () => {
    const state = createMoleculeState({
      elements: ['H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 2, y: 0, z: 0 },
      ],
    });
    expect(state.cell?.enabled).toBe(false);
    state.toggleCellVisibilityEnhanced();
    expect(state.showCell).toBe(true);
    expect(state.cell.enabled).toBe(true);
    expect(state.cell.a.x).toBeCloseTo(2 + 2, 6); // 1 Å padding on each side
    state.positions[0].x = -999;
    state.markPositionsChanged();
    const minX = state.cell.originOffset.x - 1e-6;
    const maxX = state.cell.originOffset.x + state.cell.a.x + 1e-6;
    expect(state.positions[0].x).toBeGreaterThanOrEqual(minX);
    expect(state.positions[0].x).toBeLessThan(maxX);
  });
});
