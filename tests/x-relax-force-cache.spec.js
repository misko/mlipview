/** @jest-environment jsdom */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {}, stopRenderLoop: () => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: { add: () => {} },
    },
    camera: { attachControl: () => {} },
  }),
}));

describe('x-relax-force-cache', () => {
  test.skip('relaxStep seeds cache and geometry change triggers new compute (TODO migrate to pure WS stub)', async () => {
    jest.resetModules();
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      getState: jest.fn(() => ({ connected: true })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => {}),
      requestSingleStep: jest.fn(async () => ({
        positions: [[0, 0, 0]],
        forces: [[0.05, 0, 0]],
        energy: -4.9,
      })),
      waitForEnergy: jest.fn(async () => ({
        energy: -4.9,
        forces: [[0.05, 0, 0]],
      })),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setTestHook: jest.fn(),
      ack: jest.fn(),
    };
    jest.unstable_mockModule('../public/fairchem_ws_client.js', () => ({
      getWS: () => wsStub,
    }));
    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.createElement('canvas');
    canvas.id = 'viewer';
    canvas.getBoundingClientRect = () => ({ left: 0, top: 0 });
    document.body.appendChild(canvas);
    const api = await initNewViewer(canvas, { elements: ['O'], positions: [[0, 0, 0]], bonds: [] });

    await api.relaxStep();
    const version0 = api.getForceCacheVersion();

    await api.ff.computeForces({ sync: true });
    expect(api.getForceCacheVersion()).toBe(version0);

    api.state.positions[0].x += 0.1;
    api.state.markPositionsChanged();

    await api.ff.computeForces({ sync: true });
    expect(api.getForceCacheVersion()).toBeGreaterThan(version0);
  });
});
