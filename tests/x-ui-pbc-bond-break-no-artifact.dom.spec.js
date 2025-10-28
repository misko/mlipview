/** @jest-environment jsdom */

import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

function createManipulation(state) {
  return {
    rotateBond({ j, theta }) {
      const p = state.positions[j];
      state.positions[j] = { x: p.x + Math.abs(theta) * 2.5, y: p.y, z: p.z };
      state.markPositionsChanged();
      return true;
    },
  };
}

describe('x-ui-pbc-bond-break-no-artifact', () => {
  test('breaking bond under PBC clears ghost instances and hides masters', () => {
    const cell = {
      a: { x: 10, y: 0, z: 0 },
      b: { x: 0, y: 10, z: 0 },
      c: { x: 0, y: 0, z: 10 },
      originOffset: { x: 0, y: 0, z: 0 },
      enabled: true,
    };
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.4, y: 0, z: 0 },
      ],
      cell,
    });
    state.showCell = true;
    state.showGhostCells = true;

    const bondService = createBondService(state);
    bondService.recomputeAndStore();

    const view = createMoleculeView({ onPointerObservable: { add() {} } }, state);
    view.rebuildBonds();
    view.rebuildGhosts();

    const manipulation = createManipulation(state);
    manipulation.rotateBond({ i: 0, j: 1, theta: 2.0 });
    bondService.recomputeAndStore();
    view.rebuildBonds();
    view.rebuildGhosts();

    let ghostInstances = 0;
    let anyVisible = false;
    for (const group of view._internals.ghostBondGroups.values()) {
      ghostInstances += group.mats.length;
      anyVisible = anyVisible || group.master.isVisible;
    }
    for (const group of view._internals.ghostAtomGroups.values()) {
      if (group.mats.length === 0) {
        anyVisible = anyVisible || group.master.isVisible;
      }
    }
    expect(ghostInstances).toBe(0);
    expect(anyVisible).toBe(false);
  });
});
