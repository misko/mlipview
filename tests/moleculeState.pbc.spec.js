import { createMoleculeState } from '../public/domain/moleculeState.js';

// Provide minimal BABYLON mocks required by imports from test setup

describe('moleculeState PBC integration', () => {
  test('toggleCellVisibilityEnhanced synthesizes 1Å padded cell and wraps on move', () => {
    const state = createMoleculeState({ elements:['H','H'], positions:[ {x:0,y:0,z:0}, {x:2,y:0,z:0} ] });
    // Initially cell disabled
    expect(state.cell?.enabled).toBe(false);
    // Toggle PBC on -> synthetic cell with padding
    state.toggleCellVisibilityEnhanced();
    expect(state.showCell).toBe(true);
    expect(state.cell.enabled).toBe(true);
    // Length along a should be dx + 2*1
    expect(state.cell.a.x).toBeCloseTo(2 + 2*1, 6);
    // Move an atom outside and ensure markPositionsChanged wraps it back
    state.positions[0].x = -999;
    state.markPositionsChanged();
    expect(state.positions[0].x).toBeGreaterThanOrEqual(state.cell.originOffset.x - 1e-6);
    const maxX = state.cell.originOffset.x + state.cell.a.x + 1e-6;
    expect(state.positions[0].x).toBeLessThan(maxX);
  });
});
