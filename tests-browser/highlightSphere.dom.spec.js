/**
 * Basic DOM test ensuring highlight meshes start hidden when no selection.
 */
let createMoleculeState, createMoleculeView;
beforeAll(async () => {
  ({ createMoleculeState } = await import('../public/domain/moleculeState.js'));
  ({ createMoleculeView } = await import('../public/render/moleculeView.js'));
});

function makeScene() {
  return {
    onBeforeRenderObservable: {
      _c: [],
      add(fn) {
        this._c.push(fn);
      },
      run() {
        this._c.forEach((f) => f());
      },
    },
  };
}

describe('highlight sphere hidden initially', () => {
  test('no selection => both highlight meshes hidden', () => {
    const st = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.4, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const scene = makeScene();
    const view = createMoleculeView(scene, st);
    scene.onBeforeRenderObservable.run();
    const hi = view._internals.highlight;
    expect(hi.atom.isVisible).toBe(false);
    expect(hi.bond.isVisible).toBe(false);
  });
});
