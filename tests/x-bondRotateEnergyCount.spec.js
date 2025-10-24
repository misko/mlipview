/**
 * @jest-environment jsdom
 */

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {} },
    scene: {
      onPointerObservable: { add: () => {} },
      render: () => {},
    },
    camera: { attachControl: () => {} },
  }),
}));

if (!global.BABYLON) {
  global.BABYLON = {
    Vector3: class {
      constructor(x = 0, y = 0, z = 0) {
        this.x = x;
        this.y = y;
        this.z = z;
      }
    },
    Color4: class {
      constructor(r, g, b, a) {
        this.r = r;
        this.g = g;
        this.b = b;
        this.a = a;
      }
    },
    ArcRotateCamera: class {},
    HemisphericLight: class {},
    Engine: class {},
  };
}

function setupDOM() {
  const viewer = document.createElement('canvas');
  viewer.id = 'viewer';
  document.body.appendChild(viewer);
  const energyWrapper = document.createElement('div');
  energyWrapper.id = 'energyPlot';
  document.body.appendChild(energyWrapper);
  const energyCanvas = document.createElement('canvas');
  energyCanvas.id = 'energyCanvas';
  energyCanvas.width = 260;
  energyCanvas.height = 80;
  energyCanvas.getContext = () => ({
    clearRect() {},
    beginPath() {},
    moveTo() {},
    lineTo() {},
    stroke() {},
    arc() {},
    fill() {},
    fillStyle: null,
    strokeStyle: null,
  });
  energyWrapper.appendChild(energyCanvas);
  const energyLabel = document.createElement('div');
  energyLabel.id = 'energyLabel';
  energyWrapper.appendChild(energyLabel);
  return viewer;
}

describe('x-bond rotate energy count', () => {
  test('rotateBond does not add energy samples', async () => {
    const canvas = setupDOM();
    const mod = await import('../public/index.js');
    const api = await mod.initNewViewer(canvas, {
      elements: [6, 6, 1],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.4, y: 0, z: 0 },
        { x: 2.8, y: 0, z: 0 },
      ],
      bonds: [
        { i: 0, j: 1 },
        { i: 1, j: 2 },
      ],
    });

    const before = api.debugEnergySeriesLength();
    api.debugSelectBond({ i: 0, j: 1, index: 0 });
    const changed = api.manipulation.rotateBond(Math.PI / 10);
    expect(changed).toBe(true);
    await new Promise((r) => setTimeout(r, 80));
    const after = api.debugEnergySeriesLength();
    expect(after - before).toBe(0);
  });
});
