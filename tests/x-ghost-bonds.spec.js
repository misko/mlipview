import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('ghost cell toggles', () => {
  test('toggling cell and ghosts updates flags and versions', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.5, y: 0, z: 0 },
      ],
      cell: {
        a: { x: 5, y: 0, z: 0 },
        b: { x: 0, y: 5, z: 0 },
        c: { x: 0, y: 0, z: 5 },
        enabled: true,
        originOffset: { x: 0, y: 0, z: 0 },
      },
    });
    const cellVersion0 = state.versions.cell;
    state.toggleCellVisibility();
    state.toggleGhostCells();
    expect(state.showCell).toBe(true);
    expect(state.showGhostCells).toBe(true);
    expect(state.versions.cell).toBeGreaterThan(cellVersion0);
    state.toggleGhostCells();
    expect(state.showGhostCells).toBe(false);
  });
});
