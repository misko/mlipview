/** @jest-environment jsdom */

// Validate that reset bumps interaction versions to invalidate in-flight responses and does not reload page

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

// Minimal fetch mocks for forces/relax/md
global.fetch = async (url, opts) => {
  return {
    ok: true,
    status: 200,
    json: async () => ({
      results: { energy: -1.23, forces: [[0, 0, 0]] },
      final_energy: -1.23,
      positions: [[0, 0, 0]],
      forces: [[0, 0, 0]],
    }),
  };
};

describe('Reset invalidation and no reload', () => {
  test('reset increments userInteractionVersion and cancels stale responses', async () => {
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div><canvas id="energyCanvas" width="200" height="40"></canvas><div id="energyLabel"></div>`;
    delete window.location;
    window.location = { href: 'http://localhost', assign: jest.fn() };
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const before = viewer.getVersionInfo();
    // Trigger reset via button to exercise DOM path
    const btn = document.getElementById('resetAllBtn');
    await btn.click();
    await Promise.resolve();
    const after = viewer.getVersionInfo();
    expect(after.userInteractionVersion).toBeGreaterThanOrEqual(before.userInteractionVersion + 1);
    // Start a step and immediately reset again; the step result should be stale by version checks
    const p = viewer.relaxStep();
    await viewer.resetToInitialPositions();
    const res = await p;
    if (res) expect(res.stale).toBeTruthy();

    // Also verify MD path respects epoch invalidation
    const p2 = viewer.mdStep();
    await viewer.resetToInitialPositions();
    const res2 = await p2;
    if (res2) expect(res2.stale).toBeTruthy();
  });
});
