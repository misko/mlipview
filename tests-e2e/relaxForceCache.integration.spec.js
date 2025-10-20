/** @jest-environment jsdom */
// Integration: relax seeds cache; compute stable until geometry changes; live WS.

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

describe('[integration] relax + force cache (live WS)', () => {
  test('relax seeds cache; compute stable; mutate -> bump', async () => {
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
    window.__MLIPVIEW_NO_AUTO_MD = true;
    window.__MLIPVIEW_TEST_MODE = true;

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
    // Use a small diatomic so idle compute is stable and cache versioning is clear
    const api = await initNewViewer(canvas, {
      elements: [6, 6],
      positions: [
        [0, 0, 0],
        [1.4, 0, 0],
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    await new Promise((r) => setTimeout(r, 100));
    // Use streaming relax for one step to seed cache
    await api.startRelaxContinuous({ maxSteps: 1 });
    await new Promise((r) => setTimeout(r, 150));
    const v0 = api.getForceCacheVersion();
    await api.ff.computeForces({ sync: true });
    const v1 = api.getForceCacheVersion();
    expect(v1).toBe(v0);
    api.state.positions[1].x += 0.1;
    api.state.markPositionsChanged();
    await new Promise((r) => setTimeout(r, 50));
    await api.ff.computeForces({ sync: true });
    const t1 = Date.now();
    let v2 = api.getForceCacheVersion();
    while (v2 <= v1 && Date.now() - t1 < 15000) {
      await new Promise((r) => setTimeout(r, 100));
      v2 = api.getForceCacheVersion();
    }
    expect(v2).toBeGreaterThan(v1);
  }, 30000);
});
