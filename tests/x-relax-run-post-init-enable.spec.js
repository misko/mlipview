/** @jest-environment jsdom */

import https from 'https';
import { haveServer } from './helpers/server.js';

if (typeof fetch === 'undefined') {
  global.fetch = function (nodeUrl, opts = {}) {
    return new Promise((resolve, reject) => {
      try {
        const url = new URL(nodeUrl);
        const lib = url.protocol === 'https:' ? https : require('http');
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
                json: async () => JSON.parse(body || '{}'),
                text: async () => body,
              });
            });
          }
        );
        req.on('error', reject);
        if (opts.body) req.write(opts.body);
        req.end();
      } catch (err) {
        reject(err);
      }
    });
  };
}

if (!global.BABYLON) {
  global.BABYLON = {
    TransformNode: class {},
    MeshBuilder: {
      CreateCylinder: () => ({
        dispose() {},
        setEnabled() {},
        position: { set() {} },
        rotationQuaternion: null,
        scaling: {},
        isPickable: false,
        visibility: 1,
      }),
    },
    StandardMaterial: class {
      constructor() {
        this.diffuseColor = {};
        this.emissiveColor = {};
        this.specularColor = {};
      }
    },
    Color3: class {},
    Vector3: class {
      constructor(x = 0, y = 0, z = 0) {
        this.x = x;
        this.y = y;
        this.z = z;
      }
      static Up() {
        return new global.BABYLON.Vector3(0, 1, 0);
      }
    },
    Quaternion: class {
      static Identity() {
        return {};
      }
    },
    Scene: class {},
  };
}

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

async function initViewer() {
  window.__MLIPVIEW_TEST_MODE = true;
  window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  const canvas = document.createElement('canvas');
  canvas.id = 'viewer';
  document.body.appendChild(canvas);
  const energyPlot = document.createElement('div');
  energyPlot.id = 'energyPlot';
  document.body.appendChild(energyPlot);
  const energyCanvas = document.createElement('canvas');
  energyCanvas.id = 'energyCanvas';
  energyCanvas.width = 200;
  energyCanvas.height = 60;
  energyCanvas.getContext = () => ({
    clearRect() {},
    beginPath() {},
    moveTo() {},
    lineTo() {},
    stroke() {},
    arc() {},
    fill() {},
  });
  energyPlot.appendChild(energyCanvas);
  const energyLabel = document.createElement('div');
  energyLabel.id = 'energyLabel';
  energyPlot.appendChild(energyLabel);
  const { initNewViewer } = await import('../public/index.js');
  return initNewViewer(canvas, {
    elements: [{ Z: 8 }, { Z: 1 }, { Z: 1 }],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 0.95, y: 0, z: 0 },
      { x: -0.24, y: 0.93, z: 0 },
    ],
    bonds: [],
  });
}

describe('x-relax-run-post-init-enable', () => {
  test('RELAX_LOOP flag can be toggled after init', async () => {
    if (!(await haveServer())) {
      console.warn('[skip] backend not available for relax feature test');
      return;
    }
    const api = await initViewer();
    const first = await api.startRelaxContinuous({ maxSteps: 5 });
    expect(first.disabled).not.toBe(true);
    if (typeof first.steps === 'number') {
      expect(first.steps).toBeGreaterThan(0);
    }
    api.enableFeatureFlag('RELAX_LOOP', false);
    const disabled = await api.startRelaxContinuous({ maxSteps: 3 });
    expect(disabled.disabled).toBe(true);
    api.enableFeatureFlag('RELAX_LOOP', true);
    const second = await api.startRelaxContinuous({ maxSteps: 3 });
    expect(second.disabled).not.toBe(true);
    if (typeof second.steps === 'number') {
      expect(second.steps).toBeGreaterThan(0);
    }
  }, 20000);
});
