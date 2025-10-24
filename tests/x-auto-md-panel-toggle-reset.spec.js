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

describe('x-auto MD panel toggles', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
  });

  test('MD and Relax toggles default to Off', async () => {
    const wsStub = createWsStub();
    let initNewViewer;
    let buildDesktopPanel;

    await jest.isolateModulesAsync(async () => {
      const wsModule = await import('../public/fairchem_ws_client.js');
      jest.spyOn(wsModule, 'getWS').mockReturnValue(wsStub);
      ({ initNewViewer } = await import('../public/index.js'));
      ({ buildDesktopPanel } = await import('../public/ui/desktopPanel.js'));
    });

    window.__MLIPVIEW_TEST_MODE = true;
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;

    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: ['O'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;

    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const mdToggle = document.getElementById('toggleMD');
    const relaxToggle = document.getElementById('toggleRelax');
    expect(mdToggle).toBeTruthy();
    expect(relaxToggle).toBeTruthy();
    expect(mdToggle.getAttribute('data-on')).toBe('false');
    expect(relaxToggle.getAttribute('data-on')).toBe('false');

    delete window.viewerApi;
  });
});
