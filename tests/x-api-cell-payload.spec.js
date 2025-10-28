/** @jest-environment jsdom */

describe('x-api-cell-payload', () => {
  test('USER_INTERACTION payload includes cell after enabling PBC', async () => {
    jest.resetModules();

    let seq = 0;
    let lastUserInteraction = null;
    let testHook = null;

    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      setCounters: jest.fn(),
      setTestHook: jest.fn((fn) => {
        testHook = fn;
      }),
      userInteraction: jest.fn((payload) => {
        lastUserInteraction = payload;
        if (typeof testHook === 'function') {
          testHook({ userInteraction: payload, seq: ++seq });
        }
        return ++seq;
      }),
      waitForClientSeq: jest.fn(async () => true),
      waitForEnergy: jest.fn(async () => ({ energy: -1, forces: [[[0, 0, 0]]] })),
      onResult: jest.fn(() => () => {}),
      requestSingleStep: jest.fn(async () => ({
        positions: [[0, 0, 0]],
        forces: [[0, 0, 0]],
        energy: -1,
      })),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      ack: jest.fn(),
      getState: jest.fn(() => ({ connected: true })),
    };

    jest.doMock('../public/fairchem_ws_client.js', () => ({
      getWS: () => wsStub,
    }));

    jest.doMock('../public/render/scene.js', () => ({
      createScene: async () => ({
        engine: { runRenderLoop: () => {}, stopRenderLoop: () => {} },
        scene: { meshes: [], render: () => {}, onPointerObservable: { add() {} } },
        camera: { attachControl: () => {} },
      }),
    }));

    const { initNewViewer } = await import('../public/index.js');

    if (typeof global.window === 'undefined') {
      global.window = { location: { search: '', protocol: 'http:', host: '127.0.0.1:4000' } };
    }
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';

    const origWS = global.WebSocket;
    class FakeWS {
      constructor() {
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen && this.onopen();
        }, 0);
      }
      send() {}
      close() {}
    }
    global.WebSocket = FakeWS;

    document.body.innerHTML = '<canvas id="viewer"></canvas>';
    const canvas = document.getElementById('viewer');
    Object.assign(canvas, {
      width: 640,
      height: 480,
      getContext: () => ({}),
      addEventListener: () => {},
      removeEventListener: () => {},
      getBoundingClientRect: () => ({ left: 0, top: 0, width: 640, height: 480 }),
    });

    try {
      const api = await initNewViewer(canvas, {
        elements: ['H', 'O', 'H'],
        positions: [
          [0, 0, 0],
          [0.96, 0, 0],
          [-0.24, 0.93, 0],
        ],
        bonds: [],
      });

      api.state.cell = {
        a: { x: 1, y: 0, z: 0 },
        b: { x: 0, y: 1, z: 0 },
        c: { x: 0, y: 0, z: 1 },
        enabled: false,
        originOffset: { x: 0, y: 0, z: 0 },
      };
      api.state.showCell = false;

      api.state.toggleCellVisibilityEnhanced();
      await api.ff.computeForces({ sync: true });

      const calls = wsStub.userInteraction.mock.calls.map(([payload]) => payload);
      const withCell = calls.find((payload) => payload && payload.cell);
      expect(withCell).toBeTruthy();
      expect(Array.isArray(withCell.cell)).toBe(true);
      expect(withCell.cell.length).toBe(3);
      expect(withCell.cell[0].length).toBe(3);
    } finally {
      global.WebSocket = origWS;
    }
  });
});
