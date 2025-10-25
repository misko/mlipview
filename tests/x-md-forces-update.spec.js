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

describe('x-md-forces-update', () => {
  test('mdStep updates state.forces from WS frame', async () => {
    jest.resetModules();
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      getState: jest.fn(() => ({ connected: true })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => true),
      waitForEnergy: jest.fn(async () => ({ energy: -3.0, forces: [[0.1, 0, 0], [0, 0, 0]] })),
      requestSingleStep: jest.fn(async () => ({
        positions: [
          [0, 0, 0],
          [1, 0, 0],
        ],
        forces: [
          [0.3, 0, 0],
          [0, 0, 0],
        ],
        energy: -3.0,
      })),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setTestHook: jest.fn(),
      ack: jest.fn(),
      onResult: jest.fn(() => () => {}),
    };

    jest.doMock('../public/fairchem_ws_client.js', () => ({
      getWS: () => wsStub,
    }));

    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.createElement('canvas');
    canvas.id = 'viewer';
    document.body.appendChild(canvas);

    const viewer = await initNewViewer(canvas, {
      elements: ['H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });

    await viewer.mdStep();
    expect(viewer.state.forces[0][0]).toBeCloseTo(0.3, 6);
  });
});

