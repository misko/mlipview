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

describe('x-md-temperature-slider', () => {
  test('slider value propagates to mdStep request', async () => {
    jest.resetModules();
    const wsStub = {
      ensureConnected: jest.fn(async () => true),
      getState: jest.fn(() => ({ connected: true })),
      setCounters: jest.fn(),
      userInteraction: jest.fn(() => 1),
      waitForClientSeq: jest.fn(async () => true),
      waitForEnergy: jest.fn(async () => ({ energy: -5, forces: [[0, 0, 0]] })),
      requestSingleStep: jest.fn(async ({ params }) => ({
        positions: [
          [0, 0, 0],
          [1, 0, 0],
        ],
        forces: [
          [0, 0, 0],
          [0, 0, 0],
        ],
        energy: -5,
        velocities: [],
        temperature: params?.temperature ?? null,
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

    document.body.innerHTML = `
      <canvas id="viewer"></canvas>
      <div class="hud"></div>
    `;

    const { initNewViewer } = await import('../public/index.js');
    const { initTemperatureSlider } = await import('../public/ui/temperatureSlider.js');

    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: ['O', 'H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.96, y: 0, z: 0 },
        { x: -0.24, y: 0.93, z: 0 },
      ],
      bonds: [],
    });
    window.viewerApi = viewer;
    initTemperatureSlider({ hudEl: document.querySelector('.hud'), getViewer: () => viewer });

    const slider = document.getElementById('mdTempSlider');
    slider.value = slider.max;
    slider.dispatchEvent(new Event('input'));
    const target = window.__MLIP_TARGET_TEMPERATURE;

    await viewer.mdStep();
    expect(wsStub.requestSingleStep).toHaveBeenCalled();
    const call = wsStub.requestSingleStep.mock.calls.pop();
    expect(call[0].params.temperature).toBe(target);
  });
});

