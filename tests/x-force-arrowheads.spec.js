import { createMoleculeView } from '../public/render/moleculeView.js';
import { ensureBabylonStub } from './utils/babylonStub.js';

ensureBabylonStub();

function createBus() {
  const listeners = {};
  return {
    on(ev, fn) {
      (listeners[ev] = listeners[ev] || []).push(fn);
    },
    emit(ev) {
      (listeners[ev] || []).forEach((fn) => fn());
    },
  };
}

describe('x-force-arrowheads', () => {
  test('force group populates shaft/head matrices', () => {
    const scene = { meshes: [], onBeforeRenderObservable: { add: () => {} } };
    const state = {
      elements: ['O', 'H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
        { x: -1, y: 0, z: 0 },
      ],
      bonds: [],
      bus: createBus(),
      showForces: true,
    };
    state.forces = [
      [0.2, 0, 0],
      [0.1, 0, 0],
      [-0.15, 0, 0],
    ];

    const view = createMoleculeView(scene, state);
    state.bus.emit('forcesChanged');
    const group =
      view._internals.forceGroups.get('force') || (view.rebuildForces(), view._internals.forceGroups.get('force'));
    expect(group).toBeTruthy();
    expect(group.shaftMaster._buffers['matrix'].length).toBeGreaterThan(0);
    expect(group.headMaster._buffers['matrix'].length).toBeGreaterThan(0);
  });
});
