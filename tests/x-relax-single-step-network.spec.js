/** @jest-environment jsdom */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {}, stopRenderLoop: () => {} },
    scene: { render: () => {}, onPointerObservable: { add: () => {} } },
    camera: {},
  }),
}));

describe('x-relax-single-step-network', () => {
  test.skip('relaxStep performs single requestSingleStep call (TODO migrate to pure WS stub)', async () => {
    jest.resetModules();
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      getState: jest.fn(() => ({ connected: true })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => {}),
      requestSingleStep: jest.fn(async ({ type }) => {
        if (type !== 'relax') throw new Error('unexpected type ' + type);
        return {
          positions: [[0, 0, 0]],
          forces: [[0, 0, 0]],
          energy: -9.9,
        };
      }),
      waitForEnergy: jest.fn(async () => ({ energy: -9.9, forces: [[0, 0, 0]] })),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setTestHook: jest.fn(),
      ack: jest.fn(),
    };
    jest.unstable_mockModule('../public/fairchem_ws_client.js', () => ({
      getWS: () => wsStub,
    }));

    document.body.innerHTML = '<canvas id="viewer"></canvas><div class="hud"></div>';
    const { initNewViewer } = await import('../public/index.js');
    const api = await initNewViewer(document.getElementById('viewer'), {
      elements: ['O'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });

    await api.relaxStep();
    expect(wsStub.requestSingleStep).toHaveBeenCalledTimes(1);
    expect(wsStub.requestSingleStep.mock.calls[0][0]?.type).toBe('relax');
  });
});
