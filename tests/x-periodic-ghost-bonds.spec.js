/**
 * @jest-environment jsdom
 */

import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('periodic ghost bonds without simulation', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
    document.body.innerHTML = '<canvas id="viewer" width="400" height="300"></canvas>';
    if (typeof window !== 'undefined') {
      window.__MLIPVIEW_TEST_MODE = true;
      window.__MLIP_CONFIG = window.__MLIP_CONFIG || {};
    }
  });

  afterEach(() => {
    if (typeof window !== 'undefined') {
      delete window.__MLIPVIEW_TEST_MODE;
    }
  });

  test('ghost bonds appear after enabling periodic cell while paused', async () => {
    await jest.unstable_mockModule('../public/render/scene.js', () => ({
      createScene: async () => {
        const engine = {
          runRenderLoop: () => {},
          stopRenderLoop: () => {},
          getRenderingCanvas: () => document.getElementById('viewer'),
        };
        const scene = {
          meshes: [],
          onPointerObservable: { add() {} },
          onBeforeRenderObservable: { add() {} },
          render: () => {},
          dispose: () => {},
          getEngine: () => engine,
        };
        const camera = {
          attachControl: () => {},
          detachControl: () => {},
        };
        return { engine, scene, camera };
      },
    }));

    const wsStub = {
      connect: jest.fn(async () => {}),
      ensureConnected: jest.fn(async () => true),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => {}),
      requestSingleStep: jest.fn(async () => ({})),
      waitForEnergy: jest.fn(async () => ({ energy: -1, forces: [] })),
      onResult: jest.fn(() => () => {}),
      onFrame: jest.fn(() => () => {}),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setAck: jest.fn(),
      ack: jest.fn(),
    };

    let initNewViewer;
    await jest.isolateModulesAsync(async () => {
      const wsModule = await import('../public/fairchem_ws_client.js');
      jest.spyOn(wsModule, 'getWS').mockImplementation(() => wsStub);
      ({ initNewViewer } = await import('../public/index.js'));
    });

    const canvas = document.getElementById('viewer');
    const api = await initNewViewer(canvas, {
      elements: ['C', 'C'],
      positions: [
        { x: 0.1, y: 0, z: 0 },
        { x: 2.9, y: 0, z: 0 },
      ],
      bonds: [],
      cell: {
        a: { x: 3, y: 0, z: 0 },
        b: { x: 0, y: 3, z: 0 },
        c: { x: 0, y: 0, z: 3 },
        enabled: true,
        originOffset: { x: 0, y: 0, z: 0 },
      },
    });

    api.recomputeBonds();
    expect(api.debugGhostSnapshot().ghostBondCount).toBe(0);

    if (!api.state.showCell) {
      api.state.toggleCellVisibilityEnhanced?.();
    }
    if (!api.state.showGhostCells) {
      api.state.toggleGhostCells?.();
    }

    api.recomputeBonds();

    const snapshot = api.debugGhostSnapshot();
    expect(snapshot.ghostBondCount).toBeGreaterThan(0);
    expect(Array.isArray(api.state.ghostBondMeta)).toBe(true);
    expect(api.state.ghostBondMeta.length).toBeGreaterThan(0);
  });
});
