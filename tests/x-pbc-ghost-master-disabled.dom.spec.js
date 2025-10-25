import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

function ghostInstanceCount(view) {
  return Array.from(view._internals.ghostBondGroups.values()).reduce(
    (sum, grp) => sum + (grp.mats?.length || 0),
    0
  );
}

describe('x-pbc-ghost-master-disabled', () => {
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
        { x: 1.42, y: 0, z: 0 },
      ],
      bonds: [],
      cell,
    });
    state.showCell = true;
    state.showGhostCells = true;
    const bondService = createBondService(state);
    bondService.recomputeAndStore();
    const view = createMoleculeView({ onPointerObservable: { add() {} } }, state);
    view.rebuildBonds();
    view.rebuildGhosts();

    expect(ghostInstanceCount(view)).toBeGreaterThan(0);

    state.positions[1] = { x: 6, y: 0, z: 0 };
    state.markPositionsChanged();
    bondService.recomputeAndStore();
    view.rebuildBonds();
    view.rebuildGhosts();

    let totalInstances = 0;
    let anyVisible = false;
    for (const group of view._internals.ghostBondGroups.values()) {
      totalInstances += group.mats.length;
      anyVisible = anyVisible || !!group.master.isVisible;
    }
    expect(totalInstances).toBe(0);
    expect(anyVisible).toBe(false);
  });
});
