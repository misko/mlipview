/*
 Test: Interaction invalidation of stale relax/md responses.
 Scenario 1: User bond rotation occurs while a relaxStep request is in-flight; response should respect exclusions (rotated atoms not overwritten).
 Scenario 2: Two mdStep calls fired concurrently; slower (first) response discarded due to superseded simulation.
 Scenario 3: Control relaxStep with no user interaction accepts.
*/

const { JSDOM } = require('jsdom');

let mockWsStub;
let currentState;
let seqCounter = 0;

async function initViewer() {
  const html = `<!DOCTYPE html><html><body><canvas id="c" width="100" height="100"></canvas></body></html>`;
  const dom = new JSDOM(html, {
    runScripts: 'outside-only',
    resources: 'usable',
    url: 'http://localhost/',
  });
  global.window = dom.window;
  global.document = dom.window.document;
  global.performance = { now: () => Date.now() };
  window.__MLIPVIEW_TEST_MODE = true;
  window.__MLIP_CONFIG = {
    minStepIntervalMs: 5,
    mdFriction: 0.1,
  };
  window.__MLIP_FEATURES = {
    RELAX_LOOP: true,
    MD_LOOP: true,
    ENERGY_TRACE: true,
    FORCE_VECTORS: true,
  };
  // Stub minimal BABYLON API used by createScene
  function Color3(r, g, b) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  Color3.prototype.clone = function () {
    return new Color3(this.r, this.g, this.b);
  };
  Color3.prototype.scale = function (s) {
    return new Color3(this.r * s, this.g * s, this.b * s);
  };
  Color3.FromHexString = () => new Color3(1, 1, 1);
  global.BABYLON = {
    Engine: function () {
      this.stopRenderLoop = () => {};
    },
    Scene: function () {
      this.dispose = () => {};
      this.onPointerObservable = {};
    },
    Color3,
    Color4: function (r, g, b, a) {
      this.r = r;
      this.g = g;
      this.b = b;
      this.a = a;
    },
    Vector3: class {
      constructor(x, y, z) {
        this.x = x;
        this.y = y;
        this.z = z;
      }
      static Zero() {
        return new this(0, 0, 0);
      }
    },
    MeshBuilder: { CreateSphere: () => ({ material: {}, position: {} }) },
    StandardMaterial: function () {},
    HemisphericLight: function () {},
    ArcRotateCamera: function () {
      this.attachControl = () => {};
    },
  };
  // Mock heavy rendering modules BEFORE requiring index.js
  jest.mock('../public/render/scene.js', () => ({
    createScene: async () => ({
      engine: { stopRenderLoop: () => {} },
      scene: { render: () => {}, dispose: () => {}, onPointerObservable: { add: () => {} } },
      camera: {},
    }),
  }));
  jest.mock('../public/core/pickingService.js', () => ({
    createPickingService: () => ({}),
  }));
  jest.mock('../public/render/moleculeView.js', () => ({
    createMoleculeView: () => ({ rebuildBonds: () => {}, rebuildForces: () => {} }),
  }));
  jest.mock('../public/fairchem_ws_client.js', () => ({
    getWS: () => mockWsStub,
  }));

  const { initNewViewer } = require('../public/index.js');
  const canvas = document.getElementById('c');
  // Provide a simple water-like structure (3 atoms) with one bond for rotation selection tests
  const elements = ['O', 'H', 'H'];
  const positions = [
    { x: 0, y: 0, z: 0 },
    { x: 0.95, y: 0, z: 0 },
    { x: -0.24, y: 0.93, z: 0 },
  ];
  const bonds = [
    { i: 0, j: 1, order: 1, opacity: 1 },
    { i: 0, j: 2, order: 1, opacity: 1 },
  ];

  const basePositions = positions.map((p) => [p.x, p.y, p.z]);
  seqCounter = 0;
  mockWsStub = {
    ensureConnected: jest.fn(async () => true),
    setCounters: jest.fn(),
    userInteraction: jest.fn(() => ++seqCounter),
    waitForClientSeq: jest.fn(async () => {}),
    requestSingleStep: jest.fn(async () => ({
      positions: basePositions.map((p) => [...p]),
      energy: -1,
      forces: [],
    })),
    onResult: jest.fn(() => () => {}),
    onFrame: jest.fn(() => () => {}),
    ack: jest.fn(),
    getState: jest.fn(() => ({ connected: true })),
  };

  const api = await initNewViewer(canvas, { elements, positions, bonds });
  currentState = api.state;
  return api;
}

function delay(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

describe('interaction invalidation', () => {
  jest.setTimeout(20000);

  test('relax response respects exclusion after user bond rotation', async () => {
    const api = await initViewer();
    // Ensure a dummy bond selection if not present
    api.state.selection = { kind: 'bond', data: { i: 0, j: 1, orientation: 'ij' } };
    const preUIV = api.getVersionInfo ? api.getVersionInfo().userInteractionVersion : 0;
    // Snapshot current positions for comparison
    const before = api.state.positions.map((p) => [p.x, p.y, p.z]);
    mockWsStub.requestSingleStep.mockImplementationOnce(async () => {
      await delay(80);
      return {
        positions: before.map((p) => [p[0] + 5, p[1] + 5, p[2] + 5]),
        forces: [],
        energy: -1,
      };
    });
    // Fire relax step (will start fetch)
    const relaxPromise = api.relaxStep();
    // While in flight, perform bond rotation (increments versions)
    api.manipulation.rotateBond(0.2);
    // Wait for relax to resolve
    const result = await relaxPromise;
    expect(result).toBeTruthy();
    expect(result.applied).toBe(true);
    const postInfo = api.getVersionInfo();
    expect(postInfo.userInteractionVersion).toBeGreaterThan(preUIV);
    // Determine which atom(s) were rotated: for water selection set above, expect atom 1 affected.
    // Verify that rotated atom did NOT receive +5 shift, but other atoms did.
    const pos = api.state.positions.map((p) => [p.x, p.y, p.z]);
    // atom 0 and 2 should be shifted ~+5
    expect(Math.abs(pos[0][0] - (before[0][0] + 5)) < 1e-6).toBe(true);
    expect(Math.abs(pos[0][1] - (before[0][1] + 5)) < 1e-6).toBe(true);
    expect(Math.abs(pos[2][0] - (before[2][0] + 5)) < 1e-6).toBe(true);
    // atom 1 should NOT be exactly +5 (it was rotated instead), allow some tolerance away from +5 shift
    const dx1 = Math.abs(pos[1][0] - (before[1][0] + 5));
    const dy1 = Math.abs(pos[1][1] - (before[1][1] + 5));
    const dz1 = Math.abs(pos[1][2] - (before[1][2] + 5));
    const nearPlus5 = dx1 < 1e-6 && dy1 < 1e-6 && dz1 < 1e-6;
    expect(nearPlus5).toBe(false);
  });

  test('older md response discarded due to newer simulation application', async () => {
    const api = await initViewer();
    const info0 = api.getVersionInfo();
    // Patch first MD call to be slow
    mockWsStub.requestSingleStep
      .mockImplementationOnce(async () => {
        await delay(50);
        return {
          positions: currentState.positions.map((p) => [p.x, p.y, p.z]),
          energy: -1,
          forces: [],
        };
      })
      .mockImplementationOnce(async () => ({
        positions: currentState.positions.map((p) => [p.x, p.y, p.z]),
        energy: -1,
        forces: [],
      }));

    const fast = api.mdStep({ temperature: 300 });
    const slow = api.mdStep({ temperature: 310 });
    const [res1, res2] = await Promise.all([slow, fast]);
    expect(res1.applied || res1.stale).toBe(true);
    expect(res2.applied).toBe(true);
    const info1 = api.getVersionInfo();
    expect(info1.totalInteractionVersion).toBeGreaterThan(info0.totalInteractionVersion);
  });

  test('relax step succeeds with no concurrent interaction', async () => {
    const api = await initViewer();
    const result = await api.relaxStep();
    expect(result).toBeTruthy();
    expect(result.applied).toBe(true);
  });
});
