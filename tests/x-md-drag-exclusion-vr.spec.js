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

describe('x-md-drag-exclusion-vr', () => {
  test('held atom retains drag position after mdStep', async () => {
    jest.resetModules();
    let hook = null;
    const listeners = [];
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      getState: jest.fn(() => ({ connected: true })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => true),
      waitForEnergy: jest.fn(async () => ({ energy: -1, forces: [[0, 0, 0]] })),
      requestSingleStep: jest.fn(async () => ({
        positions: [
          [0.4, 0, 0],
          [10, 0, 0],
        ],
        forces: [[0, 0, 0], [0, 0, 0]],
        energy: -1.2,
      })),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setTestHook: jest.fn((fn) => {
        hook = fn;
      }),
      ack: jest.fn(),
      onResult: jest.fn((fn) => {
        listeners.push(fn);
        return () => {
          const idx = listeners.indexOf(fn);
          if (idx >= 0) listeners.splice(idx, 1);
        };
      }),
      injectTestResult(frame) {
        listeners.forEach((fn) => fn(frame));
      },
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
      bonds: [],
    });

    viewer.selection.clickAtom(0);
    const intersector = () => ({ x: 3.0, y: 0, z: 0 });
    viewer.manipulation.beginDrag(intersector, {
      planePoint: { x: 0, y: 0, z: 0 },
      planeNormal: { x: 0, y: 1, z: 0 },
      source: 'vr',
    });
    viewer.manipulation.updateDrag(intersector);

    const heldBefore = { ...viewer.state.positions[0] };
    const step = viewer.mdStep();
    hook?.({ type: 'USER_INTERACTION' });
    wsStub.injectTestResult({
      positions: [
        [0.4, 0, 0],
        [10, 0, 0],
      ],
      forces: [[0, 0, 0], [0, 0, 0]],
      energy: -1.2,
    });
    await step;
    viewer.manipulation.endDrag();

    const p0 = viewer.state.positions[0];
    const p1 = viewer.state.positions[1];
    expect(p0.x).toBeCloseTo(heldBefore.x, 6);
    expect(p1.x).toBeCloseTo(10, 6);
  });
});

