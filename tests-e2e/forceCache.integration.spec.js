/** @jest-environment jsdom */
// Integration: force cache version increments only when geometry changes; live WS idle compute.

const { haveServer, getApiBase } = require('../tests/helpers/server.js');

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

describe('[integration] force cache versioning (live WS)', () => {
  test('baseline -> compute no-change -> mutate -> compute bumps', async () => {
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
      console.warn('[integration] WebSocket API not available; skipping');
      return;
    }
    const base = getApiBase();
    const u = new URL(base);
    delete window.location;
    window.location = {
      protocol: u.protocol,
      host: u.host,
      origin: u.origin,
      search: '',
      href: u.href,
    };
    window.__MLIPVIEW_SERVER = base;

    document.body.innerHTML = '';
    const canvas = document.createElement('canvas');
    canvas.id = 'viewer';
    canvas.addEventListener = () => {};
    document.body.appendChild(canvas);
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

    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(canvas, {
      elements: [6, 6],
      positions: [
        [0, 0, 0],
        [1.4, 0, 0],
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    // Seed with one relax frame via streaming to ensure cache + WS are warm
    await api.startRelaxContinuous({ maxSteps: 1 });
    await new Promise((r) => setTimeout(r, 150));

    // Allow baseline to seed cache
    await new Promise((r) => setTimeout(r, 200));
    const v0 = api.getForceCacheVersion();

    // computeForces (no geometry change) => version stable
    await api.ff.computeForces({ sync: true });
    const v1 = api.getForceCacheVersion();
    expect(v1).toBe(v0);

    // mutate then compute => version increases
    api.state.positions[1].x += 0.2;
    api.state.markPositionsChanged();
    // Allow markPositionsChanged debounced paths to settle
    await new Promise((r) => setTimeout(r, 50));
    // Use a streaming MD step to ensure a fresh energy/forces frame arrives
    await api.startMDContinuous({ steps: 1, temperature: 1200 });
    // Wait for version to reflect mutation
    const t0 = Date.now();
    let v2 = api.getForceCacheVersion();
    while (v2 <= v1 && Date.now() - t0 < 15000) {
      await new Promise((r) => setTimeout(r, 100));
      v2 = api.getForceCacheVersion();
    }
    expect(v2).toBeGreaterThan(v1);
  }, 40000);
});
