/** @jest-environment jsdom */

// Validate the HUD forces toggle button flips viewer state and button text.

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

let origFetch;

beforeAll(() => {
  origFetch = global.fetch;
  global.fetch = async () => ({ ok: true, status: 200, json: async () => ({ results: {} }) });
});

afterAll(() => {
  global.fetch = origFetch;
});

function buildHud() {
  document.body.innerHTML = `
    <canvas id="viewer"></canvas>
    <div class="hud">
      <button id="btnRelax"></button>
      <button id="btnRelaxRun"></button>
      <button id="btnMD"></button>
      <button id="btnMDRun"></button>
      <select id="forceProviderSel"></select>
      <button id="btnCell"></button>
      <button id="btnGhosts"></button>
      <button id="btnToggleForces"></button>
      <span id="status">Ready</span>
      <select id="xrModeSelect"></select>
    </div>
    <div id="energyPlot">
      <canvas id="energyCanvas" width="260" height="80"></canvas>
      <div id="energyLabel"></div>
    </div>`;
}

describe('x-forces-toggle', () => {
  test('button toggles viewer.showForces flag', async () => {
    buildHud();
    const { initNewViewer } = await import('../public/index.js');
    const { initForcesToggle } = await import('../public/ui/forcesToggle.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }, { Z: 1 }, { Z: 1 }],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.96, y: 0, z: 0 },
        { x: -0.24, y: 0.93, z: 0 },
      ],
      bonds: [],
    });
    window.viewerApi = viewer;
    initForcesToggle({ getViewer: () => viewer });

    const btn = document.getElementById('btnToggleForces');
    expect(btn).toBeTruthy();
    expect(viewer.state.showForces).toBe(false);
    expect(btn.textContent.toLowerCase()).toContain('forces on');

    btn.click();
    expect(viewer.state.showForces).toBe(true);
    expect(btn.textContent.toLowerCase()).toContain('forces off');

    btn.click();
    expect(viewer.state.showForces).toBe(false);
    expect(btn.textContent.toLowerCase()).toContain('forces on');
  });
});
