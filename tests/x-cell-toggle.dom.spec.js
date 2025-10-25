/** @jest-environment jsdom */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => {
    const engine = {
      runRenderLoop: () => {},
      stopRenderLoop: () => {},
      getRenderingCanvas: () =>
        (globalThis.document && globalThis.document.getElementById('viewer')) || null,
    };
    const scene = {
      meshes: [],
      render: () => {},
      onBeforeRenderObservable: { add: () => {} },
      onPointerObservable: { add: () => {} },
      pointerX: 0,
      pointerY: 0,
      pick: () => ({ hit: false }),
      createPickingRay: () => ({ direction: { x: 0, y: 0, z: 1 } }),
      getEngine: () => engine,
    };
    const camera = {
      attachControl: () => {},
      detachControl: () => {},
      position: { x: 0, y: 0, z: -10 },
      inertialAlphaOffset: 0,
      inertialBetaOffset: 0,
      inertialRadiusOffset: 0,
      inertialPanningX: 0,
      inertialPanningY: 0,
    };
    return { engine, scene, camera };
  },
}));

const originalFetch = global.fetch;

beforeEach(() => {
  document.body.innerHTML =
    '<canvas id="viewer"></canvas><div class="hud"></div><span id="status">Ready</span>';
  global.fetch = async () =>
    new Response(
      JSON.stringify({
        results: { energy: -1.0, forces: [[0, 0, 0]] },
      }),
      { status: 200 }
    );
});

afterEach(() => {
  global.fetch = originalFetch;
});

describe('x-cell-toggle', () => {
  test('button toggles cell/ghosts and label text', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const { initCellToggle } = await import('../public/ui/cellToggle.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');

    const canvas = document.getElementById('viewer');
    const viewer = await initNewViewer(canvas, {
      elements: ['O', 'H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.96, y: 0, z: 0 },
        { x: -0.24, y: 0.93, z: 0 },
      ],
      bonds: [],
    });
    window.viewerApi = viewer;
    window.__MLIP_DEV_MODE = true;
    buildDesktopPanel({ attachTo: document.querySelector('.hud') });
    initCellToggle({ getViewer: () => viewer, hudEl: document.querySelector('.hud') });

    const btn = document.getElementById('btnCell');
    expect(btn).toBeTruthy();
    expect(viewer.state.showCell).toBe(false);
    expect(viewer.state.showGhostCells).toBe(false);
    expect(btn.textContent.toLowerCase()).toContain('cell on');

    btn.click();
    expect(viewer.state.showCell).toBe(true);
    expect(viewer.state.showGhostCells).toBe(true);
    expect(btn.textContent.toLowerCase()).toContain('cell off');

    btn.click();
    expect(viewer.state.showCell).toBe(false);
    expect(viewer.state.showGhostCells).toBe(false);
    expect(btn.textContent.toLowerCase()).toContain('cell on');
  });
});
