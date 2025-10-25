import { createMoleculeView } from '../public/render/moleculeView.js';

function createBus() {
  const map = new Map();
  return {
    on: (event, fn) => {
      (map.get(event) || map.set(event, []).get(event)).push(fn);
    },
    emit: (event) => {
      (map.get(event) || []).forEach((fn) => fn());
    },
  };
}

describe('x-molecule-switch-highlight', () => {
  test('switching to single atom hides bond highlight and disables empty masters', async () => {
    const molState = {
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.4, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1, opacity: 1 }],
      selection: { kind: 'bond', data: { i: 0, j: 1 } },
      showCell: false,
      showGhostCells: false,
      cell: { enabled: false },
      bus: createBus(),
    };
    const view = createMoleculeView({}, molState);
    const hi = view._internals.highlight;
    molState.bus.emit('selectionChanged');
    expect(hi.bond.isVisible).toBe(true);

    molState.elements = ['C'];
    molState.positions = [{ x: 0, y: 0, z: 0 }];
    molState.bonds = [];
    molState.bus.emit('bondsChanged');

    expect(hi.bond.isVisible).toBe(false);
    for (const [, group] of view._internals.bondGroups) {
      if (group.mats.length === 0) {
        const master = group.master;
        expect(master.isVisible).toBe(false);
      }
    }
  });
});
