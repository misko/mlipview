/** @jest-environment jsdom */

describe('x-md-initial-temperature-ws', () => {
  test('startMDContinuous propagates XYZ temperature into WS startSimulation', async () => {
    jest.resetModules();

    let seq = 0;
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => (++seq)),
      waitForClientSeq: jest.fn(async () => true),
      onResult: jest.fn(() => () => {}),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      waitForEnergy: jest.fn(async () => ({ energy: -1, forces: [] })),
      requestSingleStep: jest.fn(),
      ack: jest.fn(),
      setTestHook: jest.fn(),
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

    const { parseXYZ } = await import('../public/util/xyzLoader.js');
    const { applyParsedToViewer } = await import('../public/util/moleculeLoader.js');
    const { initNewViewer } = await import('../public/index.js');

    document.body.innerHTML = '<canvas id="viewer"></canvas><div id="hud"></div>';
    const canvas = document.getElementById('viewer');
    Object.assign(canvas, {
      width: 800,
      height: 600,
      getContext: () => ({}),
      addEventListener: () => {},
      removeEventListener: () => {},
      getBoundingClientRect: () => ({ left: 0, top: 0, width: 800, height: 600 }),
    });

    window.__MLIPVIEW_TEST_MODE = true;
    window.__MLIPVIEW_NO_AUTO_MD = true;

    const api = await initNewViewer(canvas, { elements: [], positions: [], bonds: [] });
    const xyz = `3\n temp=500K\nH 0 0 0\nH 1 0 0\nO 0 1 0\n`;
    const parsed = parseXYZ(xyz);
    applyParsedToViewer(api, parsed);

    await api.startMDContinuous({ steps: 5 });
    await Promise.resolve();

    expect(wsStub.startSimulation).toHaveBeenCalled();
    const payload = wsStub.startSimulation.mock.calls[0][0];
    expect(payload.type).toBe('md');
    expect(payload.params.temperature).toBe(500);
  });
});

