import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.ts';
import { createMoleculeView } from '../public/render/moleculeView.js';

function dist(a, b) {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  const dz = a.z - b.z;
  return Math.hypot(dx, dy, dz);
}

describe('x-pbc-long-bond-crossing-x', () => {
  test('long primaries across X have zero opacity while ghosts handle short bond', () => {
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
        { x: 0.1, y: 0.1, z: 0.1 },
        { x: 9.9, y: 0.1, z: 0.1 },
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

    let sawLong = 0;
    let visibleLong = 0;
    for (const group of view._internals.bondGroups.values()) {
      for (const idx of group.indices) {
        const L = dist(state.positions[idx.i], state.positions[idx.j]);
        if (L > 5) {
          sawLong++;
          const alpha = idx.opacity != null ? idx.opacity : 1;
          if (alpha > 1e-3) visibleLong++;
        }
      }
    }
    expect(sawLong).toBeGreaterThanOrEqual(1);
    expect(visibleLong).toBe(0);

    const ghostCount = Array.from(view._internals.ghostBondGroups.values()).reduce(
      (sum, group) => sum + group.mats.length,
      0
    );
    expect(ghostCount).toBeGreaterThan(0);
  });
});
