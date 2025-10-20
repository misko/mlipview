/*
Verify that running mdStep / relaxStep with no user interactions does NOT increment userInteractionVersion.
We allow totalInteractionVersion to increment on applied simulation results.
*/

const { JSDOM } = require('jsdom');

async function initViewer() {
  const html = `<!DOCTYPE html><html><body><canvas id="c" width="50" height="50"></canvas></body></html>`;
  const dom = new JSDOM(html, {
    runScripts: 'outside-only',
    resources: 'usable',
    url: 'http://localhost/',
  });
  global.window = dom.window;
  global.document = dom.window.document;
  global.performance = { now: () => Date.now() };
  window.__MLIPVIEW_TEST_MODE = true;
  // Lightweight mocks
  jest.mock('../public/render/scene.js', () => ({
    createScene: async () => ({
      engine: { stopRenderLoop: () => {} },
      scene: { render: () => {}, dispose: () => {}, onPointerObservable: { add: () => {} } },
      camera: {},
    }),
  }));
  jest.mock('../public/core/pickingService.js', () => ({ createPickingService: () => ({}) }));
  jest.mock('../public/render/moleculeView.js', () => ({
    createMoleculeView: () => ({ rebuildBonds: () => {}, rebuildForces: () => {} }),
  }));
  const { initNewViewer } = require('../public/index.js');
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
  return await initNewViewer(document.getElementById('c'), { elements, positions, bonds });
}

// Patch fetch to fast-return synthetic MD/Relax responses
function mockFetchMDRelax(api) {
  const orig = global.fetch || window.fetch;
  global.fetch = async (url, opts) => {
    const u = String(url);
    if (u.includes('/serve/relax')) {
      return new Response(
        JSON.stringify({
          positions: api.state.positions.map((p) => [p.x + 0.01, p.y, p.z]),
          forces: [],
          final_energy: -10,
        }),
        { status: 200, headers: { 'Content-Type': 'application/json' } }
      );
    }
    if (u.includes('/serve/md')) {
      return new Response(
        JSON.stringify({
          positions: api.state.positions.map((p) => [p.x + 0.02, p.y, p.z]),
          forces: [],
          final_energy: -11,
        }),
        { status: 200, headers: { 'Content-Type': 'application/json' } }
      );
    }
    if (u.includes('/serve/simple')) {
      return new Response(
        JSON.stringify({
          results: {
            energy: -100,
            forces: [
              [0, 0, 0],
              [0, 0, 0],
              [0, 0, 0],
            ],
          },
        }),
        { status: 200, headers: { 'Content-Type': 'application/json' } }
      );
    }
    return orig(url, opts);
  };
  return () => {
    global.fetch = orig;
  };
}

if (typeof Response === 'undefined') {
  global.Response = class {
    constructor(body, init) {
      this._body = body;
      this.status = init.status;
      this.headers = init.headers || {};
    }
    async json() {
      return JSON.parse(this._body);
    }
    async text() {
      return this._body;
    }
    get ok() {
      return this.status >= 200 && this.status < 300;
    }
  };
}

describe('userInteractionVersion isolation', () => {
  test('mdStep & relaxStep do not bump userInteractionVersion', async () => {
    const api = await initViewer();
    const cleanup = mockFetchMDRelax(api);
    // Establish baseline energy
    await api.baselineEnergy();
    const before = api.getVersionInfo();
    // Run several simulation steps without any user actions
    for (let i = 0; i < 3; i++) await api.mdStep();
    for (let i = 0; i < 2; i++) await api.relaxStep();
    // Allow debounce timer (50ms) to fire if it were going to
    await new Promise((r) => setTimeout(r, 80));
    const after = api.getVersionInfo();
    expect(after.userInteractionVersion).toBe(before.userInteractionVersion);
    expect(after.totalInteractionVersion).toBeGreaterThanOrEqual(
      before.totalInteractionVersion + 1
    );
    cleanup();
  });
});
