/** @jest-environment jsdom */

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

async function setupViewer() {
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
  await Promise.resolve();
  return viewer;
}

describe('x-selection-periodic-click', () => {
  test('clicking periodic table entry clears selection and shows static properties', async () => {
    const viewer = await setupViewer();
    viewer.selection.clickAtom(0);
    const oxygenCell = document.querySelector('#miniPeriodic .pt-el[data-symbol=\"O\"]');
    expect(oxygenCell).toBeTruthy();
    oxygenCell.dispatchEvent(new MouseEvent('click', { bubbles: true }));

    expect((document.getElementById('selElementName').textContent || '').toLowerCase()).toContain(
      'oxygen'
    );
    expect(document.getElementById('selPosition').textContent).toBe('(-,-,-)');
    expect(document.getElementById('selAtomicWeight').textContent).toBe('15.999');
    expect(document.getElementById('selSphereA').style.display).toBe('block');
    expect(document.getElementById('selSphereB').style.display).toBe('none');
    expect(viewer.selection.get().kind).toBeNull();
  });
});
