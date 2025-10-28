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
      multiPick: undefined,
      pick: () => ({ hit: false }),
      pointerX: 0,
      pointerY: 0,
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
  document.body.innerHTML = `
    <canvas id="viewer"></canvas>
    <div id="app"></div>
  `;
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

describe('x-desktop-bond-pick', () => {
  test('clicking bond mesh selects bond and rotates', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');

    window.__MLIP_DEV_MODE = true;
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });

    const bondGroups = viewer.view._internals.bondGroups;
    const master = Array.from(bondGroups.values())[0]?.master;
    expect(master).toBeTruthy();

    viewer.scene.pick = () => ({
      hit: true,
      pickedMesh: master,
      thinInstanceIndex: 0,
    });

    const evt =
      typeof PointerEvent !== 'undefined'
        ? new PointerEvent('pointerdown', { bubbles: true, clientX: 20, clientY: 20 })
        : new Event('pointerdown', { bubbles: true });
    document.getElementById('viewer').dispatchEvent(evt);

    expect(viewer.selection.get().kind).toBe('bond');
    const highlight = viewer.view._internals.highlight;
    expect(highlight?.bond?.isVisible).toBe(true);

    const rotateBtns = document.getElementById('rotateBtns');
    const plus = document.getElementById('bondRotPlus');
    expect(rotateBtns?.style.display).toBe('inline-flex');
    const version0 = viewer.state.versions.positions;
    plus.click();
    await Promise.resolve();
    expect(viewer.state.versions.positions).toBeGreaterThan(version0);
  });
});
