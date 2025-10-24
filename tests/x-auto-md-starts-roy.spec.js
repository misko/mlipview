/** @jest-environment jsdom */

import { jest } from '@jest/globals';
import fs from 'fs';
import path from 'path';
import { installBabylonStub } from './helpers/installBabylonStub.js';
import { createWsStub } from './helpers/createWsStub.js';

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {}, stopRenderLoop: () => {} },
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

describe('x-auto MD starts ROY', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
    delete window.__MLIPVIEW_TEST_MODE;
    delete window.__MLIPVIEW_NO_AUTO_MD;
  });

  afterEach(() => {
    delete window.viewerApi;
  });

  test('auto-starts MD and streams frames', async () => {
    const wsStub = createWsStub();

    let initNewViewer;
    let loadDefault;
    await jest.isolateModulesAsync(async () => {
      const wsModule = await import('../public/fairchem_ws_client.js');
      jest.spyOn(wsModule, 'getWS').mockReturnValue(wsStub);
      ({ initNewViewer } = await import('../public/index.js'));
      ({ loadDefault } = await import('../public/util/moleculeLoader.js'));
    });

    const royPath = path.resolve(process.cwd(), 'public/molecules/roy.xyz');
    const royText = fs.readFileSync(royPath, 'utf8');
    const origFetch = global.fetch;
    global.fetch = async (url) => {
      if (typeof url === 'string' && /molecules\/roy\.xyz$/.test(url)) {
        return { ok: true, status: 200, text: async () => royText };
      }
      return { ok: false, status: 404, text: async () => 'not found' };
    };

    try {
      document.body.innerHTML = '';
      const canvas = document.createElement('canvas');
      canvas.id = 'viewer';
      canvas.getContext = () => ({
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
      document.body.appendChild(canvas);

      const hud = document.createElement('div');
      hud.className = 'hud';
      document.body.appendChild(hud);
      const status = document.createElement('span');
      status.id = 'status';
      hud.appendChild(status);
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

      const viewer = await initNewViewer(canvas, { elements: [], positions: [], bonds: [] });
      window.viewerApi = viewer;

      const loaded = await loadDefault(viewer);
      expect(loaded.file).toBe('molecules/roy.xyz');

      await new Promise((resolve) => setTimeout(resolve, 50));
      expect(wsStub.startSimulation).toHaveBeenCalled();

      const N = viewer.state.positions.length || 1;
      const basePos = viewer.state.positions.map((p) => [p.x, p.y, p.z]);
      for (let i = 0; i < 20; i++) {
        const pos = basePos.map((p, idx) => [p[0] + i * 0.0001 + (idx === 0 ? i * 0.00005 : 0), p[1], p[2]]);
        const forces = Array.from({ length: N }, (_, k) => [0.001 * (((k + i) % 3) + 1), 0, 0]);
        wsStub.emit({
          positions: pos,
          forces,
          energy: -100 - i,
          temperature: 300 + i * 0.1,
        });
      }

      await new Promise((resolve) => setTimeout(resolve, 20));

      expect(viewer.debugEnergySeriesLength()).toBeGreaterThanOrEqual(20);
    } finally {
      global.fetch = origFetch;
    }
  });
});
