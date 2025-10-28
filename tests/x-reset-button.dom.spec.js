/** @jest-environment jsdom */

// Ensure the desktop reset button triggers the viewer reset pathway instead of a hard reload.

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

describe('x-reset-button', () => {
  beforeEach(() => {
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
    delete window.location;
    window.location = {
      href: 'http://localhost',
      assign: jest.fn(),
    };
  });

  test('clicking reset button calls viewer.resetToInitialPositions()', async () => {
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const btn = document.getElementById('resetAllBtn');
    expect(btn).toBeTruthy();
    const resetSpy = jest.spyOn(viewer, 'resetToInitialPositions').mockResolvedValue(true);

    btn.click();
    await Promise.resolve();

    expect(resetSpy).toHaveBeenCalledTimes(1);
    expect(window.location.assign).not.toHaveBeenCalled();
  });
});
