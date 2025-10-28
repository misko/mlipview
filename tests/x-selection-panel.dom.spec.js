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

function wait(ms) {
  return new Promise((resolve) => setTimeout(resolve, ms));
}

async function setupViewer() {
  window.__MLIPVIEW_TEST_MODE = true;
  document.body.innerHTML = `<canvas id="viewer"></canvas><div id="app"></div>`;
  const { initNewViewer } = await import('../public/index.js');
  const { buildDesktopPanel } = await import('../public/ui/desktopPanel.js');
  const viewer = await initNewViewer(document.getElementById('viewer'), {
    elements: ['N', 'O'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
    ],
    bonds: [{ i: 0, j: 1 }],
  });
  window.viewerApi = viewer;
  buildDesktopPanel({ attachTo: document.getElementById('app') });
  await wait(0);
  return viewer;
}

describe('x-selection-panel', () => {
  test('initial state renders empty selection UI', async () => {
    await setupViewer();
    const section = document.getElementById('section-selection');
    expect(section).toBeTruthy();
    expect(document.getElementById('selSphereA').style.display).toBe('none');
    expect(document.getElementById('selSphereB').style.display).toBe('none');
    expect(document.getElementById('selElementName').textContent).toBe('—');
    expect(document.getElementById('selPosition').textContent).toBe('(-,-,-)');
  });

  test('atom selection populates details and highlights periodic cell', async () => {
    const viewer = await setupViewer();
    viewer.selection.clickAtom(0);
    const elName = document.getElementById('selElementName').textContent.toLowerCase();
    expect(elName).toContain('nitrogen');
    expect(document.getElementById('selSphereA').style.display).toBe('block');
    expect(document.getElementById('selSphereB').style.display).toBe('none');
    const nCell = document.querySelector('#miniPeriodic .pt-el[data-symbol=\"N\"]');
    const oCell = document.querySelector('#miniPeriodic .pt-el[data-symbol=\"O\"]');
    expect(nCell.classList.contains('highlight')).toBe(true);
    expect(oCell.classList.contains('highlight')).toBe(false);
  });

  test('bond selection shows pair info and length', async () => {
    const viewer = await setupViewer();
    viewer.selection.clickBond({ i: 0, j: 1, index: 0, key: '0-1' });
    const name = document.getElementById('selElementName').textContent.toLowerCase();
    expect(name.includes('nitrogen') && name.includes('oxygen')).toBe(true);
    const bondLen = document.getElementById('bondLength').textContent;
    expect(bondLen).toContain('Å');
    expect(document.getElementById('selSphereA').style.display).toBe('block');
    expect(document.getElementById('selSphereB').style.display).toBe('block');
  });

  test('clearing selection hides spheres and resets labels', async () => {
    const viewer = await setupViewer();
    viewer.selection.clickAtom(0);
    viewer.selection.clear();
    expect(document.getElementById('selSphereA').style.display).toBe('none');
    expect(document.getElementById('selSphereB').style.display).toBe('none');
    expect(document.getElementById('selElementName').textContent).toBe('—');
    const highlighted = document.querySelectorAll('#miniPeriodic .pt-el.highlight');
    expect(highlighted.length).toBe(0);
  });
});
