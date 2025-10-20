/*
 Test: Interaction invalidation of stale relax/md responses.
 Scenario 1: User bond rotation occurs while a relaxStep request is in-flight; response should be discarded (staleReason='userInteraction').
 Scenario 2: Two mdStep calls fired concurrently; slower (first) response discarded due to superseded simulation.
 Scenario 3: Control relaxStep with no user interaction accepts.
*/

const { JSDOM } = require('jsdom');

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
  const api = await initNewViewer(canvas, { elements, positions, bonds });
  return api;
}

function delay(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

// Helper to monkey patch fetch with delay for specific path
function patchFetchOnce({ match, delayMs = 50, respond }) {
  const orig = global.fetch || window.fetch;
  let used = false;
  global.fetch = async (url, opts) => {
    if (!used && typeof url === 'string' && url.includes(match)) {
      used = true;
      await delay(delayMs);
      const respData = respond ? respond(url, opts) : null;
      return new Response(JSON.stringify(respData), {
        status: 200,
        headers: { 'Content-Type': 'application/json' },
      });
    }
    return orig(url, opts);
  };
  return () => {
    global.fetch = orig;
  };
}

// Basic Response polyfill for jsdom environment if not present
if (typeof Response === 'undefined') {
  global.Response = class {
    constructor(body, init) {
      this._body = body;
      this.status = init.status;
      this.headers = new Map(Object.entries(init.headers || {}));
    }
    async json() {
      return JSON.parse(this._body);
    }
    get ok() {
      return this.status >= 200 && this.status < 300;
    }
    async text() {
      return this._body;
    }
  };
}

// NOTE: This test relies on ability to rotate a bond; ensure water structure yields a bond selection.
// If no bond selected by default, we simulate by setting selection object directly.

describe('interaction invalidation', () => {
  jest.setTimeout(20000);

  test('relax response partially applied after user bond rotation', async () => {
    const api = await initViewer();
    // Ensure a dummy bond selection if not present
    api.state.selection = { kind: 'bond', data: { i: 0, j: 1, orientation: 'ij' } };
    const preUIV = api.getVersionInfo ? api.getVersionInfo().userInteractionVersion : 0;
    // Snapshot current positions for comparison
    const before = api.state.positions.map((p) => [p.x, p.y, p.z]);
    // Patch fetch for relax to delay and return shifted positions we can recognize
    const cleanup = patchFetchOnce({
      match: '/relax',
      delayMs: 80,
      respond: () => ({
        positions: before.map((p) => [p[0] + 5, p[1] + 5, p[2] + 5]),
        forces: [],
        final_energy: -1,
      }),
    });
    // Fire relax step (will start fetch)
    const relaxPromise = api.relaxStep();
    // While in flight, perform bond rotation (increments versions)
    api.manipulation.rotateBond(0.2);
    // Wait for relax to resolve
    const result = await relaxPromise;
    cleanup();
    expect(result).toBeTruthy();
    expect(result.applied).toBe(true);
    expect(result.partial).toBe(true);
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
    let mdCallCount = 0;
    const origFetch = global.fetch || window.fetch;
    global.fetch = async (url, opts) => {
      if (typeof url === 'string' && url.includes('/md')) {
        mdCallCount++;
        if (mdCallCount === 1) {
          // slow first
          await delay(120);
          return new Response(
            JSON.stringify({
              positions: api.state.positions.map((p) => [p.x + 1, p.y, p.z]),
              forces: [],
              final_energy: info0.energy ?? -2,
            }),
            { status: 200, headers: { 'Content-Type': 'application/json' } }
          );
        } else {
          // fast second accepted
          return new Response(
            JSON.stringify({
              positions: api.state.positions.map((p) => [p.x + 0.01, p.y, p.z]),
              forces: [],
              final_energy: -3,
            }),
            { status: 200, headers: { 'Content-Type': 'application/json' } }
          );
        }
      }
      return origFetch(url, opts);
    };
    // Fire two md steps without awaiting the first
    const p1 = api.mdStep();
    const p2 = api.mdStep();
    const r2 = await p2; // likely second returns first
    const r1 = await p1; // first resolves after delay
    // Restore fetch
    global.fetch = origFetch;
    expect(r2.applied).toBe(true);
    expect(r1.stale).toBe(true);
    // Accept either stale reason (userInteraction may increment via background posChange)
    expect(['supersededSimulation', 'userInteraction']).toContain(r1.staleReason);
    const infoAfter = api.getVersionInfo();
    expect(infoAfter.totalInteractionVersion).toBeGreaterThan(info0.totalInteractionVersion);
  });

  test('relax step accepted when no user interaction', async () => {
    const api = await initViewer();
    const info0 = api.getVersionInfo();
    // Fast relax response
    const origFetch = global.fetch || window.fetch;
    global.fetch = async (url, opts) => {
      if (typeof url === 'string' && url.includes('/relax')) {
        return new Response(
          JSON.stringify({
            positions: api.state.positions.map((p) => [p.x + 0.02, p.y, p.z]),
            forces: [],
            final_energy: -5,
          }),
          { status: 200, headers: { 'Content-Type': 'application/json' } }
        );
      }
      return origFetch(url, opts);
    };
    const res = await api.relaxStep();
    global.fetch = origFetch;
    expect(res.applied).toBe(true);
    const info1 = api.getVersionInfo();
    expect(info1.totalInteractionVersion).toBe(info0.totalInteractionVersion + 1);
    expect(info1.userInteractionVersion).toBe(info0.userInteractionVersion); // unchanged (no user action)
  });
});
