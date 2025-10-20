// Test switching molecules clears bond highlight and hides empty bond masters.

describe('molecule switch clears highlight and hides empty masters', () => {
  test('switch from two-atom to single-atom', async () => {
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
    const mod = await import('../public/render/moleculeView.js');
    const view = mod.createMoleculeView({}, molState);
    const hi = view._internals.highlight;
    // Simulate selection event so highlight appears
    molState.bus.emit('selectionChanged');
    expect(hi.bond.isVisible).toBe(true);
    // Now switch to single atom molecule (no bonds)
    molState.elements = ['C'];
    molState.positions = [{ x: 0, y: 0, z: 0 }];
    molState.bonds = [];
    molState.bus.emit('bondsChanged');
    // Highlight should be hidden
    expect(hi.bond.isVisible).toBe(false);
    // All bond group masters (if any) with zero instances should be disabled/hidden
    for (const [, g] of view._internals.bondGroups) {
      if (g.mats.length === 0) {
        expect(g.master.isVisible === false || (g.master.isEnabled && !g.master.isEnabled())).toBe(
          true
        );
      }
    }
  });
});

function createBus() {
  const map = new Map();
  return {
    on: (e, f) => {
      (map.get(e) || map.set(e, []).get(e)).push(f);
    },
    emit: (e) => {
      (map.get(e) || []).forEach((f) => f());
    },
  };
}
