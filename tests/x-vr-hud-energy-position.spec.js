import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

if (!global.window) global.window = { location: { search: '' } };
if (!global.document) {
  global.document = {
    createElement: (tag) =>
      tag === 'canvas'
        ? { width: 0, height: 0, getContext: () => ({}), toDataURL: () => '' }
        : {},
  };
}

if (!global.BABYLON) global.BABYLON = {};

const B = global.BABYLON;
const vec = (x = 0, y = 0, z = 0) => ({
  x,
  y,
  z,
  add(o) {
    return vec(x + o.x, y + o.y, z + o.z);
  },
  subtract(o) {
    return vec(x - o.x, y - o.y, z - o.z);
  },
  scale(s) {
    return vec(x * s, y * s, z * s);
  },
  normalize() {
    const L = Math.hypot(x, y, z) || 1;
    return vec(x / L, y / L, z / L);
  },
  copyFrom(o) {
    this.x = o.x;
    this.y = o.y;
    this.z = o.z;
    return this;
  },
});

if (!B.Axis) B.Axis = { Z: vec(0, 0, 1), Y: vec(0, 1, 0) };
if (!B.Mesh) B.Mesh = { BILLBOARDMODE_NONE: 0 };
B.MeshBuilder = B.MeshBuilder || {};
B.MeshBuilder.CreatePlane = (name, opts = {}, scene) => {
  const plane = {
    name,
    opts,
    scene,
    position: vec(),
    rotation: { y: 0 },
    scaling: { x: 1, y: 1, z: 1 },
    billboardMode: 0,
    isPickable: true,
    metadata: {},
    dispose() {},
  };
  plane.position.copyFrom = (other) => {
    plane.position.x = other.x;
    plane.position.y = other.y;
    plane.position.z = other.z;
    return plane.position;
  };
  return plane;
};

if (!B.GUI) B.GUI = {};
B.GUI.Control = B.GUI.Control || { HORIZONTAL_ALIGNMENT_CENTER: 0 };
B.GUI.AdvancedDynamicTexture = class {
  constructor() {
    this._rootContainer = { children: [] };
  }
  static CreateForMesh(mesh, width, height) {
    const tex = new B.GUI.AdvancedDynamicTexture();
    tex._linkedMesh = mesh;
    tex.width = width;
    tex.height = height;
    return tex;
  }
  addControl(ctrl) {
    this._rootContainer.children.push(ctrl);
  }
};

B.GUI.StackPanel = class {
  constructor() {
    this.children = [];
    this.isVertical = true;
    this.horizontalAlignment = 0;
  }
  addControl(ctrl) {
    this.children.push(ctrl);
  }
};

B.GUI.Rectangle = class {
  constructor() {
    this.children = [];
    this.background = '';
    this.thickness = 0;
  }
  addControl(ctrl) {
    this.children.push(ctrl);
  }
};

B.GUI.TextBlock = class {
  constructor(id, text) {
    this.id = id;
    this.text = text;
    this.color = '#fff';
    this.fontSize = 16;
  }
};

B.GUI.Image = class {
  constructor(id, src) {
    this.id = id;
    this.src = src;
    this.width = '';
    this.height = '';
    this.stretch = 0;
  }
  markAsDirty() {}
};
B.GUI.Image.STRETCH_UNIFORM = 0;

B.GUI.Button = {
  CreateSimpleButton(id, text) {
    return {
      id,
      text,
      width: '',
      height: '',
      background: '',
      onPointerDownObservable: { add() {} },
      onPointerUpObservable: { add() {} },
      onPointerEnterObservable: { add() {} },
      onPointerOutObservable: { add() {} },
    };
  },
};

function makeScene() {
  return {
    activeCamera: {
      fov: Math.PI / 2,
      position: vec(),
      getDirection(axis) {
        if (axis === B.Axis.Z) return vec(0, 0, 1);
        if (axis === B.Axis.Y) return vec(0, 1, 0);
        return vec();
      },
    },
    onBeforeRenderObservable: {
      add() {},
    },
    getEngine() {
      return {
        getCaps() {
          return { maxTextureSize: 2048 };
        },
      };
    },
  };
}

describe('x-vr HUD energy position', () => {
  test('distance multiplier is honored', () => {
    const scene = makeScene();
    const result = ensureWorldHUD({ scene, topEnergyDistanceMult: 2.75 });
    expect(result.topEnergyDistanceMult).toBeCloseTo(2.75, 5);
  });
});
