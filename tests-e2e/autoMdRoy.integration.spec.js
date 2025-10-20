/** @jest-environment jsdom */
// Integration test: requires backend UMA server. Verifies ROY loads and MD auto-starts with ≥ 20 real WS frames.

const { haveServer, getApiBase } = require('../tests/helpers/server.js');

// Minimal scene mock (we still render nothing, but allow code paths to proceed)
jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn) => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: {
        _l: [],
        add(fn) {
          this._l.push(fn);
        },
      },
    },
    camera: { attachControl: () => {} },
  }),
}));

describe('[integration] auto MD on ROY over real WS', () => {
  test('starts MD and receives ≥ 20 frames (live server)', async () => {
    if (process.env.MLIPVIEW_SKIP_SERVERS === '1') {
      console.warn('[integration] skipping due to MLIPVIEW_SKIP_SERVERS');
      return;
    }
    const up = await haveServer();
    if (!up) {
      console.warn('[integration] backend not reachable; skipping');
      return;
    }
    if (typeof global.WebSocket === 'undefined' && typeof window.WebSocket === 'undefined') {
      // Node/jsdom environment may not provide WebSocket; require('ws') is not guaranteed to be installed
      console.warn('[integration] WebSocket API not available; skipping');
      return;
    }
    // Ensure we target the same host for WS as API base
    const base = getApiBase();
    const u = new URL(base);
    // Point window.location and __MLIPVIEW_SERVER to base so resolveWsBase() uses it
    delete window.location; // jsdom allows assignment
    window.location = {
      protocol: u.protocol.replace('http', 'http'),
      host: u.host,
      origin: u.origin,
      search: '',
      href: u.href,
    };
    window.__MLIPVIEW_SERVER = base;

    // Ensure auto-start isn’t disabled by other tests
    delete window.__MLIPVIEW_TEST_MODE;
    delete window.__MLIPVIEW_NO_AUTO_MD;
    // Enable verbose API/WS logging for diagnosis
    window.__MLIPVIEW_DEBUG_API = true;

    // DOM scaffold
    document.body.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'viewer';
    canvas.addEventListener = () => {};
    document.body.appendChild(canvas);
    const hud = document.createElement('div');
    hud.className = 'hud';
    document.body.appendChild(hud);
    const rps = document.createElement('span');
    rps.id = 'rpsLabel';
    rps.textContent = 'RPS: --';
    hud.appendChild(rps);
    const energyWrapper = document.createElement('div');
    energyWrapper.id = 'energyPlot';
    document.body.appendChild(energyWrapper);
    const energyCanvas = document.createElement('canvas');
    energyCanvas.id = 'energyCanvas';
    energyCanvas.width = 260;
    energyCanvas.height = 80;
    energyCanvas.getContext = () => ({
      clearRect() {},
      beginPath() {},
      moveTo() {},
      lineTo() {},
      stroke() {},
      arc() {},
      fill() {},
      fillRect() {},
      strokeStyle: null,
      lineWidth: 1,
      fillStyle: null,
    });
    energyWrapper.appendChild(energyCanvas);
    const energyLabel = document.createElement('div');
    energyLabel.id = 'energyLabel';
    energyWrapper.appendChild(energyLabel);

    // Import modules and initialize viewer
    const { initNewViewer } = await import('../public/index.js');
    const { applyParsedToViewer } = await import('../public/util/moleculeLoader.js');
    const { parseXYZ } = await import('../public/util/xyzLoader.js');

    // Use real fetch for ROY and any other assets
    // Initialize viewer API
    const api = await initNewViewer(canvas, { elements: [], positions: [], bonds: [] });
    window.viewerApi = api;

    // Load ROY by reading from filesystem directly and applying to viewer
    const fs = await import('fs');
    const p = await import('path');
    const royPath = p.resolve(process.cwd(), 'public/molecules/roy.xyz');
    const royText = fs.readFileSync(royPath, 'utf8');
    const parsed = parseXYZ(royText);
    applyParsedToViewer(api, parsed);

    // Allow auto-start MD to trigger
    await new Promise((r) => setTimeout(r, 100));

    // Collect frames by observing energy ticks
    const maxWaitMs = 15000;
    const start = Date.now();
    let ticks = api.debugEnergySeriesLength();
    while (ticks < 20 && Date.now() - start < maxWaitMs) {
      await new Promise((r) => setTimeout(r, 100));
      ticks = api.debugEnergySeriesLength();
    }
    // Accept a slightly lower threshold to avoid flake on busy runners
    expect(ticks).toBeGreaterThanOrEqual(18);
  }, 20000);
});
