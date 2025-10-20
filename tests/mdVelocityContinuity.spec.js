// Verifies velocity continuity and precomputed coexistence for MD steps.
import { jest } from '@jest/globals';
import { getWS } from '../public/fairchem_ws_client.js';

const API_BASE = 'http://localhost:8000';

// Minimal mocks similar to existing tests (avoid full Babylon dependency)
if (!global.window) global.window = {};
var window = global.window;
// Provide minimal location for WS base URL resolution
if (!window.location)
  window.location = { protocol: 'http:', host: 'localhost:8000', origin: 'http://localhost:8000' };
if (!global.document) global.document = { getElementById: () => null };
if (!global.BABYLON) global.BABYLON = {};
// Minimal Engine/Scene + camera/vector mocks used by render/scene.js
class MockEngine {
  constructor() {
    this._loaders = [];
  }
  runRenderLoop(fn) {
    /* noop */
  }
  stopRenderLoop() {}
  dispose() {}
}
class MockScene {
  constructor() {
    this.onBeforeRenderObservable = { add: () => {} };
    this.onPointerObservable = { add: () => {} };
    this.meshes = [];
    this._engine = new MockEngine();
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
global.BABYLON.Engine = MockEngine;
global.BABYLON.Scene = MockScene;
global.BABYLON.Vector3 = Vector3;
global.BABYLON.Color4 = Color4;
global.BABYLON.ArcRotateCamera = ArcRotateCamera;

// Provide stable elementToZ mapping indirectly via index.js import.

// Sequence storage
const requests = [];
let mdCall = 0;

// Deterministic velocities for test
function makeVel(step) {
  // 3 atoms: pattern increments slightly each step
  return [
    [0.1 + 0.01 * step, 0.0, 0.0],
    [0.0, 0.2 + 0.01 * step, 0.0],
    [0.0, 0.0, 0.3 + 0.01 * step],
  ];
}

// No fetch; simulate MD frames via WS.

// Lightweight engine/scene mocks used by index.js internals
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

describe('md velocity continuity', () => {
  test('second MD step sends velocities from first response and includes precomputed', async () => {
    const { initNewViewer } = await import('../public/index.js');
    // Stub WS
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
    // Inject a frame right after each START_SIMULATION request is sent by mdStep
    let simCount = 0;
    ws.setTestHook((msg) => {
      try {
        if (msg && msg.simulationParams) {
          const step = simCount++;
          // Deliver a frame with positions, forces, energy, velocities
          const positions = api?.state?.positions
            ? api.state.positions.map((p) => [p.x, p.y, p.z])
            : [
                [0, 0, 0],
                [0.95, 0, 0],
                [-0.24, 0.93, 0],
              ];
          const velocities = makeVel(step);
          const forces = [
            [0, 0, 0],
            [0, 0, 0],
            [0, 0, 0],
          ];
          const energy = -10.0 - 0.05 * step;
          // Ensure listener for requestSingleStep is already attached (it is before START_SIMULATION send)
          ws.injectTestResult({ positions, velocities, forces, energy, temperature: 300 });
        }
      } catch {}
    });
    const api = await initNewViewer(
      { addEventListener: () => {}, getBoundingClientRect: () => ({ left: 0, top: 0 }) },
      {
        elements: ['O', 'H', 'H'],
        positions: [
          [0, 0, 0],
          [0.95, 0, 0],
          [-0.24, 0.93, 0],
        ],
        bonds: [],
      }
    );

    // First MD step -> hook injects velocities v0
    await api.mdStep();
    // Second MD step -> hook injects velocities v1
    await api.mdStep();

    // Confirm internal state updated to second response velocities
    const stateV = api.state.dynamics.velocities;
    const secondResponseV = makeVel(1);
    expect(stateV.length).toBe(secondResponseV.length);
    for (let i = 0; i < stateV.length; i++) {
      for (let k = 0; k < 3; k++) expect(stateV[i][k]).toBeCloseTo(secondResponseV[i][k], 6);
    }
    ws.setTestHook(null);
    global.WebSocket = origWS;
    api.shutdown && api.shutdown();
  }, 15000);
});
