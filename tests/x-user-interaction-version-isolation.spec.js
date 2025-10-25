/** @jest-environment jsdom */

describe('x-user-interaction-version-isolation', () => {
  test('mdStep/relaxStep keep userInteractionVersion unchanged (WS)', async () => {
    jest.resetModules();

    const singleStepPayload = {
      positions: [[0.01, 0, 0]],
      forces: [[0, 0, 0]],
      energy: -10,
    };

    let connected = false;
    let lastSeq = 0;
    const wsStub = {
      ensureConnected: jest.fn(async () => {
        connected = true;
        return true;
      }),
      getState: jest.fn(() => ({ connected })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => {
        lastSeq += 1;
        return lastSeq;
      }),
      waitForClientSeq: jest.fn(async (seq) => {
        // mimic network delay without complex sequencing
        await Promise.resolve();
        if (!connected) throw new Error('not connected');
        if (typeof seq === 'number' && seq > lastSeq) {
          throw new Error('seq not acknowledged');
        }
        return true;
      }),
      requestSingleStep: jest.fn(async ({ type }) => {
        await Promise.resolve();
        return {
          positions: [[type === 'md' ? 0.02 : 0.03, 0, 0]],
          forces: [[0, 0, 0]],
          velocities: [[0, 0, 0]],
          temperature: type === 'md' ? 310 : undefined,
          energy: type === 'md' ? -11 : -12,
        };
      }),
      waitForEnergy: jest.fn(async () => singleStepPayload),
      startSimulation: jest.fn(),
      stopSimulation: jest.fn(),
      setTestHook: jest.fn(),
      ack: jest.fn(),
    };

    jest.doMock('../public/fairchem_ws_client.js', () => ({
      getWS: () => wsStub,
    }));

    jest.doMock('../public/render/scene.js', () => ({
      createScene: async () => ({
        engine: { runRenderLoop: () => {}, stopRenderLoop: () => {} },
        scene: { meshes: [], render: () => {}, onPointerObservable: { add() {} } },
        camera: {},
      }),
    }));

    const { initNewViewer } = await import('../public/index.js');
    document.body.innerHTML = '<canvas id="viewer"></canvas>';
    const api = await initNewViewer(document.getElementById('viewer'), {
      elements: ['O'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });

    await api.baselineEnergy();
    const before = api.getVersionInfo();

    for (let i = 0; i < 3; i++) await api.mdStep();
    for (let i = 0; i < 2; i++) await api.relaxStep();

    await new Promise((resolve) => setTimeout(resolve, 50));

    const after = api.getVersionInfo();
    expect(after.userInteractionVersion).toBe(before.userInteractionVersion);
    expect(after.totalInteractionVersion).toBeGreaterThanOrEqual(
      before.totalInteractionVersion + 1
    );
  });
});
