/** @jest-environment jsdom */

import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';
import { createWsStub } from './helpers/createWsStub.js';

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
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

describe('x-auto MD run reset', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
  });

  afterEach(() => {
    delete window.__MLIPVIEW_TEST_MODE;
    delete window.__MLIPVIEW_NO_AUTO_MD;
    delete window.viewerApi;
  });

  test('auto-started run flips button back to run', async () => {
    const startStates = [];
    const wsStub = createWsStub();
    wsStub.startSimulation = jest.fn(() => {
      const btn = document.getElementById('btnMDRun');
      if (btn) startStates.push(btn.textContent);
    });

    let initNewViewer;
    await jest.isolateModulesAsync(async () => {
      const wsModule = await import('../public/fairchem_ws_client.js');
      jest.spyOn(wsModule, 'getWS').mockReturnValue(wsStub);
      ({ initNewViewer } = await import('../public/index.js'));
    });

    document.body.innerHTML = `
      <canvas id="viewer"></canvas>
      <div class="hud">
        <button id="btnMD"></button>
        <button id="btnMDRun">run</button>
        <button id="btnRelax"></button>
        <button id="btnRelaxRun"></button>
        <select id="forceProviderSel"></select>
        <button id="btnCell"></button>
        <button id="btnGhosts"></button>
        <button id="btnToggleForces"></button>
        <span id="status">Ready</span>
        <select id="xrModeSelect"></select>
      </div>
      <div id="energyPlot"><canvas id="energyCanvas" width="260" height="80"></canvas><div id="energyLabel"></div></div>
    `;

    const energyCanvas = document.getElementById('energyCanvas');
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

    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;

    await new Promise((resolve) => setTimeout(resolve, 50));

    expect(wsStub.startSimulation).toHaveBeenCalled();
    expect(startStates).toContain('stop');

    await new Promise((resolve) => setTimeout(resolve, 0));

    expect(document.getElementById('btnMDRun').textContent).toBe('run');
  });
});
