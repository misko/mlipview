/** @jest-environment jsdom */

// Verify that a single Cell button toggles both cell and ghost states and updates label

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

beforeAll(() => {
  if (!global.BABYLON) {
    global.BABYLON = {
      TransformNode: class {},
      MeshBuilder: { CreateCylinder: () => ({}), CreateSphere: () => ({}) },
      StandardMaterial: class {},
      Color3: class {},
      Vector3: class {
        constructor(x = 0, y = 0, z = 0) {
          this.x = x;
          this.y = y;
          this.z = z;
        }
      },
      Quaternion: class {},
      Scene: class {},
      Matrix: class {},
    };
  }
});

// Minimal fetch stubs for endpoints invoked during viewer init
global.fetch = async function (url) {
  return { ok: true, status: 200, json: async () => ({}) };
};

async function setup() {
  document.body.innerHTML = `
    <canvas id="viewer"></canvas>
    <div class="hud">
      <button id="btnRelax"></button>
      <button id="btnRelaxRun"></button>
      <button id="btnMD"></button>
      <button id="btnMDRun"></button>
      <select id="forceProviderSel"></select>
      <button id="btnCell"></button>
      <button id="btnToggleForces"></button>
      <span id="status">Ready</span>
      <select id="xrModeSelect"></select>
    </div>
    <div id="energyPlot"><canvas id="energyCanvas" width="260" height="80"></canvas><div id="energyLabel"></div></div>`;
  const { initNewViewer } = await import('../public/index.js');
  const { initCellToggle } = await import('../public/ui/cellToggle.js');
  const canvas = document.getElementById('viewer');
  const viewer = await initNewViewer(canvas, {
    elements: [{ Z: 8 }, { Z: 1 }, { Z: 1 }],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 0.96, y: 0, z: 0 },
      { x: -0.24, y: 0.93, z: 0 },
    ],
    bonds: [],
  });
  window.viewerApi = viewer;
  initCellToggle({ getViewer: () => viewer });
  return viewer;
}

describe('Cell toggle (unified)', () => {
  test('single button toggles cell and ghosts and label updates', async () => {
    const v = await setup();
    const btn = document.getElementById('btnCell');
    expect(btn).toBeTruthy();
    // initial state should be OFF for both
    expect(!!v.state.showCell).toBe(false);
    expect(!!v.state.showGhostCells).toBe(false);
    expect(btn.textContent.toLowerCase()).toContain('cell on');
    // click to turn ON
    btn.click();
    expect(v.state.showCell).toBe(true);
    expect(v.state.showGhostCells).toBe(true);
    expect(btn.textContent.toLowerCase()).toContain('cell off');
    // click to turn OFF
    btn.click();
    expect(v.state.showCell).toBe(false);
    expect(v.state.showGhostCells).toBe(false);
    expect(btn.textContent.toLowerCase()).toContain('cell on');
  });
});
