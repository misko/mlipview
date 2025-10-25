/** @jest-environment jsdom */

const { haveServer, getApiBase } = require('../tests/helpers/server.js');

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: { _l: [], add(fn) { this._l.push(fn); } },
    },
    camera: { attachControl: () => {} },
  }),
}));

describe('x-md-velocity-continuity (integration)', () => {
  test('two mdStep calls populate dynamics.velocities array', async () => {
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
    const api = await initNewViewer(canvas, {
      elements: [8, 1, 1],
      positions: [
        [0, 0, 0],
        [0.95, 0, 0],
        [-0.24, 0.93, 0],
      ],
      bonds: [
        { i: 0, j: 1 },
        { i: 0, j: 2 },
      ],
    });
    await new Promise((r) => setTimeout(r, 150));

    await api.startMDContinuous({ steps: 2 });
    const t0 = Date.now();
    let velocities = api.state?.dynamics?.velocities || [];
    while ((!Array.isArray(velocities) || velocities.length !== 3) && Date.now() - t0 < 20000) {
      await new Promise((r) => setTimeout(r, 100));
      velocities = api.state?.dynamics?.velocities || [];
    }
    expect(Array.isArray(velocities)).toBe(true);
    expect(velocities.length).toBe(3);
    for (const row of velocities) {
      expect(Array.isArray(row)).toBe(true);
      expect(row.length).toBe(3);
      for (const component of row) expect(Number.isFinite(component)).toBe(true);
    }
  }, 40000);
});

