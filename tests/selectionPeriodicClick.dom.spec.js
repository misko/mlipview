/** @jest-environment jsdom */

// Verify clicking periodic table cells clears selection and shows element info with (-,-,-) position.

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

function wait(ms) {
  return new Promise((r) => setTimeout(r, ms));
}

describe('Periodic table click behavior', () => {
  async function setup() {
    window.__MLIPVIEW_TEST_MODE = true;
    document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
    const { initNewViewer } = await import('../public/index.js');
    const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
    const viewer = await initNewViewer(document.getElementById('viewer'), {
      elements: ['C'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    window.viewerApi = viewer;
    buildDesktopPanel({ attachTo: document.getElementById('app') });
    await wait(0);
    return viewer;
  }

  test('clicking O shows Oxygen and (-,-,-), clears selection', async () => {
    const v = await setup();
    // Select an atom first to ensure click clears it
    v.selection.clickAtom(0);
    const oCell = document.querySelector('#miniPeriodic .pt-el[data-symbol="O"]');
    expect(oCell).toBeTruthy();
    oCell.dispatchEvent(new MouseEvent('click', { bubbles: true }));
    const name = (document.getElementById('selElementName').textContent || '').toLowerCase();
    const pos = document.getElementById('selPosition').textContent || '';
    const w = document.getElementById('selAtomicWeight').textContent || '';
    const vdw = document.getElementById('selVdw').textContent || '';
    const a = document.getElementById('selSphereA');
    const b = document.getElementById('selSphereB');
    expect(name).toContain('oxygen');
    expect(pos).toBe('(-,-,-)');
    expect(w).toBe('15.999');
    expect(vdw).toBe('1.52');
    expect(a.style.display).toBe('block');
    expect(b.style.display).toBe('none');
  });
});
