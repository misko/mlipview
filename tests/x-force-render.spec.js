import { createMoleculeView } from '../public/render/moleculeView.js';
import { ensureBabylonStub } from './utils/babylonStub.js';

ensureBabylonStub();

function bus() {
  const map = {};
  return {
    on(ev, fn) {
      (map[ev] = map[ev] || []).push(fn);
    },
    emit(ev) {
      (map[ev] || []).forEach((fn) => fn());
    },
  };
}

describe('x-force-render', () => {
  test('force color buffer is populated with red emissive', () => {
    const scene = { meshes: [], onBeforeRenderObservable: { add: () => {} } };
    const state = {
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [],
      bus: bus(),
    };
    state.forces = [
      [0.3, 0, 0],
      [-0.2, 0, 0],
    ];
    const view = createMoleculeView(scene, state);
    state.bus.emit('forcesChanged');
    const group =
      view._internals.forceGroups.get('force') || (view.rebuildForces(), view._internals.forceGroups.get('force'));
    expect(group).toBeTruthy();
    const colors = group.master._buffers['color'];
    expect(colors).toBeTruthy();
    expect(colors[0]).toBeCloseTo(0.95, 2);
    expect(colors[1]).toBeCloseTo(0.05, 2);
    expect(colors[2]).toBeCloseTo(0.05, 2);
  });
});
