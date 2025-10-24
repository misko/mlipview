/**
 * @jest-environment jsdom
 */

import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('x-acetic acid ghosts + API payload', () => {
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

  function aceticXYZ() {
    return `32
cell: a=13.31, b=4.09, c=5.77, alpha=90.00, beta=107.0, gamma=90.00, spacegroup=P21/c, temp=400K
C 0 0 0
H 0 0 1
H 0 1 0
H 1 0 0
C 2 0 0
H 2 0 1
H 2 1 0
H 3 0 0
C 4 0 0
O 4 0 1
O 4 1 0
C 5 0 0
H 5 0 1
H 5 1 0
H 6 0 0
H 6 0 1
C 7 0 0
O 7 0 1
O 7 1 0
C 8 0 0
H 8 0 1
H 8 1 0
H 9 0 0
H 9 0 1
C 10 0 0
O 10 0 1
O 10 1 0
C 11 0 0
H 11 0 1
H 11 1 0
H 12 0 0
H 12 0 1
`;
  }

  test('ghost instances render and WS payload contains cell+pbc', async () => {
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

    const listeners = new Set();
    const wsStub = {
      connect: jest.fn(async () => {}),
      ensureConnected: jest.fn(async () => true),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => {}),
      requestSingleStep: jest.fn(async () => {
        const response = {
          positions: [],
          forces: [],
          energy: -1.0,
        };
        listeners.forEach((fn) => fn(response));
        return response;
      }),
      waitForEnergy: jest.fn(async () => ({
        energy: -1.0,
        forces: [],
      })),
      onResult: jest.fn((fn) => {
        listeners.add(fn);
        return () => listeners.delete(fn);
      }),
      onFrame: jest.fn(() => () => {}),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setAck: jest.fn(),
      ack: jest.fn(),
    };

    let initNewViewer;
    let applyParsedToViewer;
    let parseXYZ;
    await jest.isolateModulesAsync(async () => {
      const wsModule = await import('../public/fairchem_ws_client.js');
      jest.spyOn(wsModule, 'getWS').mockImplementation(() => wsStub);
      ({ initNewViewer } = await import('../public/index.js'));
      ({ applyParsedToViewer } = await import('../public/util/moleculeLoader.js'));
      ({ parseXYZ } = await import('../public/util/xyzLoader.js'));
    });

    const canvas = document.getElementById('viewer');
    const api = await initNewViewer(canvas, {
      elements: ['H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [],
    });

    const parsed = parseXYZ(aceticXYZ());
    applyParsedToViewer(api, parsed);
    if (api.state.cell) {
      api.state.cell.enabled = true;
      api.state.showCell = true;
      api.state.showGhostCells = true;
      api.state.markCellChanged?.();
    }

    await api.requestSimpleCalculateNow();

    const calls = wsStub.userInteraction.mock.calls;
    expect(calls.length).toBeGreaterThan(0);
    const payloads = calls.map((args) => args && args[0]).filter(Boolean);
    const payload = payloads.find((arg) => arg && arg.cell);
    expect(payload).toBeDefined();
    expect(payload.natoms).toBe(32);
    expect(Array.isArray(payload.cell)).toBe(true);
    expect(payload.cell.length).toBe(3);
    expect(Array.isArray(payload.cell[0])).toBe(true);

    const ghostSnapshot = api.debugGhostSnapshot();
    expect(ghostSnapshot.ghostBondCount).toBeGreaterThan(0);
  });
});
