import fs from 'fs';
import path from 'path';
import { createCalloutLayer } from '../public/render/calloutLayer.js';

function loadSn2Fixture() {
  const sn2Path = path.resolve(process.cwd(), 'public/examples/sn2/sn2.json');
  const raw = fs.readFileSync(sn2Path, 'utf8');
  return JSON.parse(raw);
}

function makeStubBabylon(plane) {
  return {
    MeshBuilder: {
      CreatePlane: () => plane,
    },
    Mesh: {
      BILLBOARDMODE_ALL: 7,
    },
    GUI: {
      AdvancedDynamicTexture: {
        CreateForMesh() {
          return {
            wrap: false,
            rootContainer: { children: [] },
            addControl(control) {
              this.rootContainer.children.push(control);
            },
            dispose() {},
          };
        },
      },
      Rectangle: class {
        constructor() {
          this.children = [];
          this.width = 1;
          this.height = 1;
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
          this.textHorizontalAlignment = 0;
          this.textVerticalAlignment = 0;
          this.color = '';
          this.textWrapping = false;
          this.fontFamily = '';
          this.fontSize = 0;
        }
      },
      Control: {
        HORIZONTAL_ALIGNMENT_CENTER: 0,
        VERTICAL_ALIGNMENT_CENTER: 0,
      },
    },
  };
}

describe('SN2 library callout offsets', () => {
  test('callout anchor.offset affects billboard position', () => {
    const sn2 = loadSn2Fixture();
    const message = sn2.timeline.controlMessages.find((m) => m.id === 'approach-phase');
    expect(message).toBeTruthy();
    const callout = message.actions.find((action) => action.type === 'overlay.callout');
    expect(callout).toBeTruthy();

    const plane = {
      billboardMode: null,
      isPickable: false,
      renderingGroupId: null,
      scaling: { x: 1, y: 1, z: 1 },
      position: { x: 0, y: 0, z: 0 },
      setEnabled: jest.fn(),
    };
    const babylon = makeStubBabylon(plane);

    const positions = sn2.viewer.positions.map(([x, y, z]) => ({ x, y, z }));
    const getAtomPositionWorld = (idx) => positions[idx] || null;
    const getBondMidpointWorld = (i, j) => {
      const a = getAtomPositionWorld(i);
      const b = getAtomPositionWorld(j);
      if (!a || !b) return null;
      return {
        x: (a.x + b.x) / 2,
        y: (a.y + b.y) / 2,
        z: (a.z + b.z) / 2,
      };
    };

    const layer = createCalloutLayer({
      scene: {},
      babylon,
      getAtomPositionWorld,
      getBondMidpointWorld,
    });

    layer.show(callout);
    const basePosition = { ...plane.position };

    const shifted = JSON.parse(JSON.stringify(callout));
    shifted.anchor.offset[1] += 3; // raise callout by 3 Å
    layer.show(shifted);

    expect(plane.position.y - basePosition.y).toBeCloseTo(3, 6);
    expect(plane.position.x).toBeCloseTo(basePosition.x, 6);
    expect(plane.position.z).toBeCloseTo(basePosition.z, 6);
  });
});
