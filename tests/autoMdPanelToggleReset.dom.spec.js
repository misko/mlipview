/** @jest-environment jsdom */

// Minimal sanity check: panel renders MD/Relax toggles and they start Off.

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

describe('Auto MD Panel Toggle Reset', () => {
  test('initial toggles are Off', async () => {
    window.__MLIPVIEW_TEST_MODE = true;
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: [{ Z: 8 }],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    const md = document.getElementById('toggleMD');
    const rx = document.getElementById('toggleRelax');
    expect(md).toBeTruthy();
    expect(rx).toBeTruthy();
    expect(md.getAttribute('data-on')).toBe('false');
    expect(rx.getAttribute('data-on')).toBe('false');
  });
});
