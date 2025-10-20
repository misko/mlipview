/** @jest-environment jsdom */
// Bond recompute visibility test:
// 1. Initialize viewer with water (O + 2 H) close enough to form 2 bonds.
// 2. Manually move hydrogens far away to break bonds; trigger markPositionsChanged + recompute.
// 3. Assert bond count goes to 0 (or < initial).
// 4. Move hydrogens back near oxygen; recompute; assert bonds restored to >=2.
// Uses mocked scene render layer.

import http from 'http';
import { haveServer } from './helpers/server.js';
import https from 'https';

if (typeof fetch === 'undefined') {
  global.fetch = function (nodeUrl, opts = {}) {
    return new Promise((resolve, reject) => {
      try {
        const url = new URL(nodeUrl);
        const lib = url.protocol === 'https:' ? https : http;
        const req = lib.request(
          url,
          { method: opts.method || 'GET', headers: opts.headers || {} },
          (res) => {
            const chunks = [];
            res.on('data', (d) => chunks.push(d));
            res.on('end', () => {
              const body = Buffer.concat(chunks).toString('utf8');
              resolve({
                ok: res.statusCode >= 200 && res.statusCode < 300,
                status: res.statusCode,
                json: async () => JSON.parse(body),
                text: async () => body,
              });
            });
          }
        );
        req.on('error', reject);
        if (opts.body) req.write(opts.body);
        req.end();
      } catch (e) {
        reject(e);
      }
    });
  };
}

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

async function setupViewer() {
  window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  // Ensure loops disabled to avoid background bond noise
  window.__MLIP_FEATURES = {
    RELAX_LOOP: false,
    MD_LOOP: false,
    ENERGY_TRACE: false,
    FORCE_VECTORS: false,
  };
  const canvas = document.createElement('canvas');
  canvas.id = 'viewer';
  document.body.appendChild(canvas);
  const energyWrapper = document.createElement('div');
  energyWrapper.id = 'energyPlot';
  document.body.appendChild(energyWrapper);
  const energyCanvas = document.createElement('canvas');
  energyCanvas.id = 'energyCanvas';
  energyCanvas.width = 300;
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
  const mod = await import('../public/index.js');
  return await mod.initNewViewer(canvas, {
    elements: [{ Z: 8 }, { Z: 1 }, { Z: 1 }],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 0.95, y: 0, z: 0 },
      { x: -0.24, y: 0.93, z: 0 },
    ],
    bonds: [],
  });
}

function distance(a, b) {
  const dx = a.x - b.x,
    dy = a.y - b.y,
    dz = a.z - b.z;
  return Math.sqrt(dx * dx + dy * dy + dz * dz);
}

describe('bond recompute on manual position changes', () => {
  test('bonds disappear when atoms separate and return when restored', async () => {
    if (!(await haveServer())) {
      console.warn('[bond recompute] server not reachable; skipping');
      return;
    }
    // Warmup simple endpoint to ensure forcefield/model loaded (mitigate first-call latency)
    try {
      await fetch('http://127.0.0.1:8000/serve/simple', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
          atomic_numbers: [8, 1, 1],
          coordinates: [
            [0, 0, 0],
            [0.95, 0, 0],
            [-0.24, 0.93, 0],
          ],
          properties: ['energy', 'forces'],
          calculator: 'uma',
        }),
      });
    } catch {}
    const api = await setupViewer();
    // Initial recompute to ensure bonds exist
    const initialBonds = api.state.bonds ? api.state.bonds.length : 0;
    expect(initialBonds).toBeGreaterThanOrEqual(1); // water should have at least one O-H bond recognized

    // Move hydrogens far away
    api.state.positions[1].x = 50;
    api.state.positions[1].y = 50;
    api.state.positions[1].z = 50;
    api.state.positions[2].x = -50;
    api.state.positions[2].y = -50;
    api.state.positions[2].z = -50;
    api.state.markPositionsChanged();
    const gone = api.bondService.recomputeAndStore();
    expect(gone.length).toBe(0);

    // Restore hydrogens near oxygen
    api.state.positions[1].x = 0.95;
    api.state.positions[1].y = 0;
    api.state.positions[1].z = 0;
    api.state.positions[2].x = -0.24;
    api.state.positions[2].y = 0.93;
    api.state.positions[2].z = 0;
    api.state.markPositionsChanged();
    const back = api.bondService.recomputeAndStore();
    expect(back.length).toBeGreaterThanOrEqual(initialBonds);

    // Also verify viewerApi.state.bonds is updated
    expect(api.state.bonds.length).toBe(back.length);
  }, 20000);
});
