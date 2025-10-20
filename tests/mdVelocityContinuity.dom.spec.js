/** @jest-environment jsdom */
// jsdom variant: validates velocity continuity using WS single-step frames
import { jest } from '@jest/globals';
import { getWS } from '../public/fairchem_ws_client.js';

function makeVel(step) {
  return [
    [0.1 + 0.01 * step, 0, 0],
    [0, 0.2 + 0.01 * step, 0],
    [0, 0, 0.3 + 0.01 * step],
  ];
}

// Provide BABYLON minimal mocks for scene creation under jsdom
if (!global.BABYLON) global.BABYLON = {};
class MockEngine {
  constructor() {}
  runRenderLoop(fn) {
    this._loop = fn;
  }
  stopRenderLoop() {}
  dispose() {}
}
class MockScene {
  constructor() {
    this.onBeforeRenderObservable = { add: () => {} };
    this.onPointerObservable = { add: () => {} };
    this.meshes = [];
  }
  render() {}
  dispose() {}
}
class Vector3 {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
}
class Color4 {
  constructor(r, g, b, a) {
    this.r = r;
    this.g = g;
    this.b = b;
    this.a = a;
  }
}
class ArcRotateCamera {
  constructor() {}
  attachControl() {}
}
Object.assign(global.BABYLON, {
  Engine: MockEngine,
  Scene: MockScene,
  Vector3,
  Color4,
  ArcRotateCamera,
  PointerEventTypes: { POINTERDOWN: 1 },
});

// Mocks
jest.unstable_mockModule('../public/render/moleculeView.js', () => ({
  createMoleculeView: () => ({ rebuildBonds: () => {}, _internals: { forceGroups: new Map() } }),
}));
jest.unstable_mockModule('../public/domain/moleculeState.js', () => ({
  createMoleculeState: ({ elements, positions }) => ({
    elements,
    positions: positions.map((p) => ({ x: p[0], y: p[1], z: p[2] })),
    bonds: [],
    bus: { on: () => {}, emit: () => {} },
    markPositionsChanged() {},
    dynamics: {},
  }),
}));
jest.unstable_mockModule('../public/domain/bondService.js', () => ({
  createBondService: () => ({ recomputeAndStore: () => [] }),
}));
jest.unstable_mockModule('../public/domain/selectionService.js', () => ({
  createSelectionService: () => ({}),
}));
jest.unstable_mockModule('../public/core/pickingService.js', () => ({
  createPickingService: () => ({}),
}));
jest.unstable_mockModule('../public/domain/manipulationService.js', () => ({
  createManipulationService: () => ({
    beginDrag: () => {},
    updateDrag: () => {},
    endDrag: () => {},
    setDragPlane: () => {},
    rotateBond: () => {},
  }),
}));
jest.unstable_mockModule('../public/vr/setup.js', () => ({
  createVRSupport: () => ({ init: async () => ({ supported: false }) }),
}));
jest.unstable_mockModule('../public/vr/vr-picker.js', () => ({ createVRPicker: () => ({}) }));
jest.unstable_mockModule('../public/fairchem_provider.js', () => ({
  createFairChemForcefield: () => ({}),
}));
jest.unstable_mockModule('../public/util/funcCount.js', () => ({ __count: () => {} }));

// No fetch; use WS client injection

describe('md velocity continuity (jsdom)', () => {
  test('velocity reuse between md steps', async () => {
    const canvas = document.createElement('canvas');
    canvas.getBoundingClientRect = () => ({ left: 0, top: 0 });
    // Stub WS and capture
    const origWS = global.WebSocket;
    class FakeWS {
      constructor() {
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen && this.onopen();
        }, 0);
      }
      send() {}
      close() {}
      onopen() {}
      onmessage() {}
      onerror() {}
    }
    global.WebSocket = FakeWS;
    const ws = getWS();
    ws.setTestHook(() => {});
    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(canvas, {
      elements: ['O', 'H', 'H'],
      positions: [
        [0, 0, 0],
        [0.95, 0, 0],
        [-0.24, 0.93, 0],
      ],
      bonds: [],
    });
    // First MD step: inject velocities v0
    const p1 = api.mdStep();
    ws.injectTestResult({
      positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
      velocities: makeVel(0),
      forces: [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
      ],
      energy: -5.0,
      temperature: 300,
    });
    await p1;
    // Second MD step: inject velocities v1
    const p2 = api.mdStep();
    ws.injectTestResult({
      positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
      velocities: makeVel(1),
      forces: [
        [0, 0, 0],
        [0, 0, 0],
        [0, 0, 0],
      ],
      energy: -5.1,
      temperature: 300,
    });
    await p2;
    // Validate state velocities continuity
    const stateV1 = api.state.dynamics.velocities;
    for (let i = 0; i < stateV1.length; i++)
      for (let k = 0; k < 3; k++) expect(stateV1[i][k]).toBeCloseTo(makeVel(1)[i][k], 6);
    global.WebSocket = origWS;
  });
});
