import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

const cell = {
  a: { x: 10, y: 0, z: 0 },
  b: { x: 0, y: 10, z: 0 },
  c: { x: 0, y: 0, z: 10 },
  originOffset: { x: 0, y: 0, z: 0 },
  enabled: true,
};

function dist(a, b) {
  const dx = a.x - b.x;
  const dy = a.y - b.y;
  const dz = a.z - b.z;
  return Math.hypot(dx, dy, dz);
}

describe('x-pbc-long-bond-crossing', () => {
  test('no long opaque primary bond spans the periodic cell', () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0.1, y: 0.1, z: 0.1 },
        { x: 9.9, y: 9.9, z: 9.9 },
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

    const longs = [];
    for (const [key, group] of view._internals.bondGroups.entries()) {
      for (const b of group.indices) {
        const L = dist(state.positions[b.i], state.positions[b.j]);
        if (L > 6) {
          longs.push({ key, i: b.i, j: b.j, len: L, opacity: b.opacity ?? 1 });
        }
      }
    }
    expect(longs.length).toBe(0);
  });
});
