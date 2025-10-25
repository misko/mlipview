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

describe('x-energy-reset', () => {
  test('resetToInitialPositions clears energy trace', async () => {
    jest.resetModules();
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      getState: jest.fn(() => ({ connected: true })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => true),
      waitForEnergy: jest.fn(async () => ({
        energy: -2.0,
        forces: [[0, 0, 0]],
      })),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      requestSingleStep: jest.fn(async () => ({
        positions: [[0, 0, 0]],
        forces: [[0, 0, 0]],
        energy: -2.0,
      })),
      setTestHook: jest.fn(),
      ack: jest.fn(),
    };

    jest.doMock('../public/fairchem_ws_client.js', () => ({
      getWS: () => wsStub,
    }));

    document.body.innerHTML = `
      <canvas id="viewer"></canvas>
      <div id="hud"></div>
      <div id="energyPlot"></div>
      <canvas id="energyCanvas" width="260" height="80"></canvas>
      <div id="energyLabel"></div>
    `;

    const { initNewViewer } = await import('../public/index.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [[0, 0, 0]],
      bonds: [],
    });

    await viewer.ff.computeForces({ sync: true });
    expect(viewer.debugEnergySeriesLength()).toBeGreaterThanOrEqual(1);

    await viewer.resetToInitialPositions();
    await Promise.resolve();
    expect(viewer.debugEnergySeriesLength()).toBe(0);
  });
});

