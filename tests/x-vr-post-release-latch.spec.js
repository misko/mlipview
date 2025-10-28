/**
 * @jest-environment jsdom
 */

import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('x-vr post-release latch', () => {
  beforeEach(() => {
    jest.resetModules();
    jest.useFakeTimers();
    installBabylonStub();
    document.body.innerHTML = '<canvas id="viewer" width="320" height="240"></canvas>';
    if (typeof window !== 'undefined') {
      window.__MLIPVIEW_TEST_MODE = true;
      window.__MLIP_CONFIG = window.__MLIP_CONFIG || {};
    }
  });

  afterEach(() => {
    jest.useRealTimers();
    if (typeof window !== 'undefined') {
      delete window.__MLIPVIEW_TEST_MODE;
    }
  });

  test('MD frame arriving after VR drag respects latch window', async () => {
    const stepQueue = [];

    await jest.unstable_mockModule('../public/render/scene.js', () => ({
      createScene: async () => {
        const engine = {
          runRenderLoop: () => {},
          stopRenderLoop: () => {},
          getRenderingCanvas: () => document.getElementById('viewer'),
        };
        const scene = {
          onPointerObservable: { add: () => {} },
          onBeforeRenderObservable: { add: () => {} },
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

    let initNewViewer;
    let queueResponse;
    const listeners = new Set();
    const wsStub = {
      connect: jest.fn(async () => {}),
      ensureConnected: jest.fn(async () => true),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 0),
      waitForClientSeq: jest.fn(async () => {}),
      requestSingleStep: jest.fn(async () => {
        if (!stepQueue.length) throw new Error('No step response queued');
        const next = stepQueue.shift();
        listeners.forEach((fn) => fn(next));
        return next;
      }),
      waitForEnergy: jest.fn(async () => ({
        energy: -1.0,
        forces: [
          [0, 0, 0],
          [0, 0, 0],
        ],
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

    await jest.isolateModulesAsync(async () => {
      const wsModule = await import('../public/fairchem_ws_client.js');
      jest.spyOn(wsModule, 'getWS').mockImplementation(() => wsStub);
      queueResponse = (resp) => stepQueue.push(resp);
      ({ initNewViewer } = await import('../public/index.js'));
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

    api.selection.clickAtom(0);
    const intersector = () => ({ x: 2.5, y: 0, z: 0 });
    expect(api.manipulation.beginDrag(intersector)).toBe(true);
    expect(api.manipulation.updateDrag(intersector)).toBe(true);
    const dragged = { ...api.state.positions[0] };

    api.manipulation.endDrag();

    queueResponse({
      positions: [
        [0, 0, 0],
        [1, 0, 0],
      ],
      forces: [
        [0, 0, 0],
        [0, 0, 0],
      ],
      velocities: [
        [0, 0, 0],
        [0, 0, 0],
      ],
      energy: -1.5,
      temperature: 300,
    });
    const first = await api.mdStep({});
    expect(first?.applied).toBe(true);

    const p0 = api.state.positions[0];
    const p1 = api.state.positions[1];
    expect(p0.x).toBeCloseTo(dragged.x, 6);
    expect(p1.x).toBeCloseTo(1, 6);

    jest.advanceTimersByTime(1000);
    queueResponse({
      positions: [
        [10, 0, 0],
        [11, 0, 0],
      ],
      forces: [
        [0, 0, 0],
        [0, 0, 0],
      ],
      velocities: [
        [0, 0, 0],
        [0, 0, 0],
      ],
      energy: -0.5,
      temperature: 300,
    });
    const second = await api.mdStep({});
    expect(second?.applied).toBe(true);

    expect(api.state.positions[0].x).toBeCloseTo(10, 6);
    expect(api.state.positions[1].x).toBeCloseTo(11, 6);
  });
});
