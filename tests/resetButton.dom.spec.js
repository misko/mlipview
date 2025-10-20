/** @jest-environment jsdom */

// Ensure the reset button is injected and clicking triggers a navigation attempt

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: (fn) => {} },
    scene: {
      meshes: [],
      render: () => {},
      onPointerObservable: {
        _l: [],
        add(fn) {
          this._l.push(fn);
        },
      },
    },
    camera: { attachControl: () => {} },
  }),
}));

global.fetch = async () => ({ ok: true, status: 200, json: async () => ({ results: {} }) });

describe('Reset button', () => {
  test('exists and onclick triggers in-app reset (no reload)', async () => {
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
    const assignSpy = jest.fn();
    delete window.location;
    // Provide minimal location with assign
    window.location = { href: 'http://localhost', assign: assignSpy };
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
    await btn.click();
    // Give any microtasks a tick
    await Promise.resolve();
    expect(resetSpy).toHaveBeenCalled();
    expect(assignSpy).not.toHaveBeenCalled();
  });
});
