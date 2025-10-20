/** @jest-environment jsdom */

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

global.fetch = async () => ({
  ok: true,
  status: 200,
  json: async () => ({
    results: { energy: -1.0, forces: [[0, 0, 0]] },
    final_energy: -1.0,
    positions: [[0, 0, 0]],
    forces: [[0, 0, 0]],
  }),
});

describe('Energy plot clears on reset', () => {
  test('debugEnergySeriesLength becomes 0 after reset', async () => {
    document.body.innerHTML = `
      <canvas id="viewer"></canvas>
      <div id="app"></div>
      <canvas id="energyCanvas" width="200" height="40"></canvas>
      <div id="energyLabel"></div>
    `;
    delete window.location;
    window.location = { href: 'http://localhost' };
    const { initNewViewer } = await import('../public/index.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;
    // Force at least one energy point to be drawn by triggering a compute (baseline already triggers one path)
    await viewer.ff.computeForces({ sync: true });
    // Attempt to plot via an interaction (if not yet plotted)
    viewer.debugRecordInteraction('test');
    // May or may not be >0 depending on mocks, but after reset it must be 0
    await viewer.resetToInitialPositions();
    expect(viewer.debugEnergySeriesLength()).toBe(0);
  });
});
