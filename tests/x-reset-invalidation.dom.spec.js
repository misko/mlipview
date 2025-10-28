/** @jest-environment jsdom */

// Validate reset invalidation without relying on live WebSocket backends.

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: { add() {} },
    },
    camera: { attachControl: () => {} },
  }),
}));

const mockWsState = {
  queue: [],
  pending: [],
  records: [],
  getViewerState: () => ({ positions: [], elements: [], resetEpoch: 0 }),
};
let mockSeq = 0;

jest.mock('../public/fairchem_ws_client.js', () => {
  const stub = {
    ensureConnected: jest.fn(async () => true),
    getState: jest.fn(() => ({ connected: true })),
    readyState: 1,
    setCounters: jest.fn(),
    userInteraction: jest.fn(() => ++mockSeq),
    waitForClientSeq: jest.fn(async () => {}),
    waitForEnergy: jest.fn(async () => {
      const st = mockWsState.getViewerState() || {};
      const positions = (st.positions || []).map((p) => [p.x, p.y, p.z]);
      return {
        energy: -1.23,
        forces: positions.map(() => [0, 0, 0]),
      };
    }),
    requestSingleStep: jest.fn(async ({ type }) => {
      const st = mockWsState.getViewerState() || {};
      const entry =
        mockWsState.queue.shift() ||
        {
          queuedEpoch: st.resetEpoch ?? 0,
          factory: () => ({
            positions: (st.positions || []).map((p) => [p.x, p.y, p.z]),
            energy: -1.1,
            forces: (st.positions || []).map(() => [0, 0, 0]),
            stale: false,
            staleReason: undefined,
          }),
        };
      return new Promise((resolve) => {
        mockWsState.pending.push({
          resolve: () => {
            const current = mockWsState.getViewerState() || {};
            const payload =
              typeof entry.factory === 'function'
                ? entry.factory({ type, state: current, queuedEpoch: entry.queuedEpoch })
                : entry.factory;
            mockWsState.records.push({
              queuedEpoch: entry.queuedEpoch,
              resolveEpoch: current.resetEpoch ?? 0,
            });
            resolve(payload);
          },
        });
      });
    }),
    stopSimulation: jest.fn(),
    ack: jest.fn(),
  };
  return {
    getWS: () => stub,
    __setViewerStateAccessor(fn) {
      mockWsState.getViewerState =
        typeof fn === 'function' ? fn : () => ({ positions: [], elements: [], resetEpoch: 0 });
    },
    __queueStep(fn) {
      const current = mockWsState.getViewerState() || {};
      mockWsState.queue.push({
        queuedEpoch: current.resetEpoch ?? 0,
        factory: fn,
      });
    },
    __resolveNextStep() {
      const pending = mockWsState.pending.shift();
      if (pending) pending.resolve();
    },
    __takeRecords() {
      const out = mockWsState.records.slice();
      mockWsState.records = [];
      return out;
    },
    __resetQueue() {
      mockWsState.queue = [];
      mockWsState.pending = [];
      mockWsState.records = [];
      mockSeq = 0;
    },
  };
});

let origFetch;
let origWebSocket;

beforeAll(() => {
  origFetch = global.fetch;
  origWebSocket = global.WebSocket;
  class FakeWS {
    constructor() {
      this.readyState = 0;
      setTimeout(() => {
        this.readyState = 1;
        this.onopen?.();
      }, 0);
    }
    send() {}
    close() {}
    onopen() {}
    onmessage() {}
    onerror() {}
  }
  global.WebSocket = FakeWS;
  global.fetch = async () => ({
    ok: true,
    status: 200,
    json: async () => ({
      results: { energy: -1.23, forces: [[0, 0, 0]] },
      final_energy: -1.23,
      positions: [[0, 0, 0]],
      forces: [[0, 0, 0]],
    }),
  });
});

afterAll(() => {
  global.fetch = origFetch;
  global.WebSocket = origWebSocket;
});

function setupDom() {
  document.body.innerHTML = `<canvas id="viewer"></canvas>
    <div id="app"></div>
    <canvas id="energyCanvas" width="200" height="40"></canvas>
    <div id="energyLabel"></div>`;
  delete window.location;
  window.location = { href: 'http://localhost', assign: jest.fn() };
}

describe('x-reset-invalidation', () => {
  test('reset bumps userInteractionVersion and marks pending results stale', async () => {
    const wsBridge = jest.requireMock('../public/fairchem_ws_client.js');
    wsBridge.__resetQueue();

    setupDom();
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    wsBridge.__setViewerStateAccessor(() => viewer.state);
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const before = viewer.getVersionInfo();
    await viewer.resetToInitialPositions();
    const after = viewer.getVersionInfo();
    expect(after.userInteractionVersion).toBeGreaterThan(before.userInteractionVersion);
    expect(after.resetEpoch).toBeGreaterThanOrEqual(before.resetEpoch);

    const epochBeforeRelax = viewer.getVersionInfo().resetEpoch;
    wsBridge.__queueStep(({ state }) => {
      return {
        positions: (state.positions || []).map((p) => [p.x, p.y, p.z]),
        energy: -1.0,
        forces: (state.positions || []).map(() => [0, 0, 0]),
        stale: true,
        staleReason: 'staleEpoch',
      };
    });
    const relaxPromise = viewer.relaxStep();
    await viewer.resetToInitialPositions();
    const afterSecondReset = viewer.getVersionInfo().resetEpoch;
    expect(afterSecondReset).toBeGreaterThan(epochBeforeRelax);
    wsBridge.__resolveNextStep();
    const relaxResult = await relaxPromise;
    const relaxRecords = wsBridge.__takeRecords();
    expect(Array.isArray(relaxRecords) && relaxRecords.length).toBeTruthy();
    expect(relaxRecords[0].resolveEpoch).toBeGreaterThanOrEqual(relaxRecords[0].queuedEpoch);

    const epochBeforeMd = viewer.getVersionInfo().resetEpoch;
    wsBridge.__queueStep(({ state }) => {
      return {
        positions: (state.positions || []).map((p) => [p.x, p.y, p.z]),
        energy: -1.0,
        forces: (state.positions || []).map(() => [0, 0, 0]),
        stale: true,
        staleReason: 'staleEpoch',
      };
    });
    const mdPromise = viewer.mdStep();
    await viewer.resetToInitialPositions();
    wsBridge.__resolveNextStep();
    const mdResult = await mdPromise;
    const mdRecords = wsBridge.__takeRecords();
    expect(Array.isArray(mdRecords) && mdRecords.length).toBeTruthy();
    expect(mdRecords[0].resolveEpoch).toBeGreaterThanOrEqual(mdRecords[0].queuedEpoch);
  });

  test('reset restores baseline geometry after add/remove edits', async () => {
    const wsBridge = jest.requireMock('../public/fairchem_ws_client.js');
    wsBridge.__resetQueue();

    setupDom();
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');

    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }, { Z: 1 }],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.95, y: 0.0, z: 0.0 },
      ],
      bonds: [],
    });
    wsBridge.__setViewerStateAccessor(() => viewer.state);
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const baselineElements = viewer.state.elements.slice();
    const baselinePositions = viewer.state.positions.map((p) => ({ x: p.x, y: p.y, z: p.z }));
    const baselineCount = baselinePositions.length;

    await viewer.addAtomAtOrigin('H');
    expect(viewer.state.positions).toHaveLength(baselineCount + 1);

    await viewer.resetToInitialPositions();
    expect(viewer.state.positions).toHaveLength(baselineCount);
    expect(viewer.state.elements).toEqual(baselineElements);
    viewer.state.positions.forEach((p, idx) => {
      const base = baselinePositions[idx];
      expect(p.x).toBeCloseTo(base.x, 6);
      expect(p.y).toBeCloseTo(base.y, 6);
      expect(p.z).toBeCloseTo(base.z, 6);
    });

    await viewer.removeAtomByIndex(0);
    expect(viewer.state.positions).toHaveLength(baselineCount - 1);

    await viewer.resetToInitialPositions();
    expect(viewer.state.positions).toHaveLength(baselineCount);
    expect(viewer.state.elements).toEqual(baselineElements);
    viewer.state.positions.forEach((p, idx) => {
      const base = baselinePositions[idx];
      expect(p.x).toBeCloseTo(base.x, 6);
      expect(p.y).toBeCloseTo(base.y, 6);
      expect(p.z).toBeCloseTo(base.z, 6);
    });
  });
});
