import { createCalloutLayer } from '../public/render/calloutLayer.js';

function createBabylonStub() {
  const createdPlanes = [];
  const babylon = {
    MeshBuilder: {
      CreatePlane: (name, opts, scene) => {
        const mesh = {
          name,
          opts,
          scene,
          scaling: { x: 1, y: 1, z: 1 },
          position: { x: 0, y: 0, z: 0 },
          setEnabled: jest.fn(),
        };
        createdPlanes.push(mesh);
        return mesh;
      },
    },
    Mesh: {
      BILLBOARDMODE_ALL: 7,
    },
    GUI: {
      AdvancedDynamicTexture: {
        CreateForMesh: () => ({
          wrap: false,
          rootContainer: { children: [] },
          addControl(control) {
            this.rootContainer.children.push(control);
          },
          dispose: jest.fn(),
        }),
      },
      Rectangle: class {
        constructor() {
          this.children = [];
          this.background = '';
          this.alpha = 1;
          this.cornerRadius = 0;
        }
        addControl(control) {
          this.children.push(control);
        }
      },
      TextBlock: class {
        constructor() {
          this.text = '';
          this.color = '';
          this.fontFamily = '';
          this.fontSize = 0;
          this.textHorizontalAlignment = 0;
          this.textVerticalAlignment = 0;
          this.textWrapping = false;
        }
      },
      Control: {
        HORIZONTAL_ALIGNMENT_CENTER: 0,
        VERTICAL_ALIGNMENT_CENTER: 0,
      },
    },
  };
  return { babylon, createdPlanes };
}

describe('calloutLayer', () => {
  const positions = [
    { x: 0, y: 0, z: 0 },
    { x: 2, y: 2, z: 0 },
    { x: -1, y: 1, z: 0 },
  ];

  const getAtomPositionWorld = (idx) => positions[idx] || null;
  const getBondMidpointWorld = (i, j) => {
    const a = getAtomPositionWorld(i);
    const b = getAtomPositionWorld(j);
    if (!a || !b) return null;
    return { x: (a.x + b.x) / 2, y: (a.y + b.y) / 2, z: (a.z + b.z) / 2 };
  };

  test('supports multiple simultaneous callouts', () => {
    const { babylon } = createBabylonStub();
    const layer = createCalloutLayer({
      scene: {},
      babylon,
      getAtomPositionWorld,
      getBondMidpointWorld,
    });

    const callouts = [
      { key: 'callout-a', text: 'A', anchor: { mode: 'atom', atoms: [0] }, offset: [0, 1, 0] },
      { key: 'callout-b', text: 'B', anchor: { mode: 'atom', atoms: [1] }, offset: [0, 2, 0] },
    ];

    layer.showAll(callouts);

    const state = layer.getState();
    expect(state.visible).toBe(true);
    expect(state.count).toBe(2);
    const texts = state.configs.map((cfg) => cfg.text).sort();
    expect(texts).toEqual(['A', 'B']);

    // Removing one callout keeps the other rendered.
    layer.showAll([callouts[1]]);
    const stateAfter = layer.getState();
    expect(stateAfter.count).toBe(1);
    expect(stateAfter.configs[0].text).toBe('B');

    layer.hide();
    expect(layer.getState().visible).toBe(false);
  });
});
