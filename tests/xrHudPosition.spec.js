/**
 * xrHudPosition.spec
 * Verifies VR/AR world HUD bar vertical positioning logic.
 *
 * Intent: The top energy plot panel should hover just above center (small positive Y offset)
 * rather than high near top; the bottom control bar should sit somewhat below center but
 * not extremely low. This test encodes the desired geometry so it FAILS with the current
 * (legacy) constants (-0.6 top, +0.6 bottom) and will PASS once we lower the top panel and
 * slightly raise the bottom panel.
 *
 * Strategy:
 *  - Provide minimal BABYLON + GUI stubs required by ensureWorldHUD.
 *  - Create a mock scene & camera (fov = PI/2, origin) and invoke ensureWorldHUD.
 *  - Trigger follow observers once (simulating a frame) to ensure final positions.
 *  - Assert top panel vertical offset is within a narrow band just above center.
 *  - Assert bottom panel vertical offset magnitude reflects a modest drop below center.
 *
 * Notes:
 *  - This is a geometry-oriented behavioral test; it does NOT aim to validate GUI rendering.
 *  - Ranges chosen allow small numeric drift but will flag large regressions.
 */

import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

// Provide a lightweight window/document environment.
if (!global.window) global.window = { location: { search: '' } };
if (!global.document) {
  global.document = {
    createElement: (tag) => {
      if (tag === 'canvas') {
        return {
          width: 0,
          height: 0,
          getContext: () => ({}),
          toDataURL: () => 'data:image/png;base64,',
        };
      }
      return { style: {} };
    },
  };
}

// Axis vectors
function makeVector(x = 0, y = 0, z = 0) {
  return {
    x,
    y,
    z,
    add(v) {
      return makeVector(this.x + v.x, this.y + v.y, this.z + v.z);
    },
    subtract(v) {
      return makeVector(this.x - v.x, this.y - v.y, this.z - v.z);
    },
    scale(s) {
      return makeVector(this.x * s, this.y * s, this.z * s);
    },
    normalize() {
      const L = Math.hypot(this.x, this.y, this.z) || 1;
      return makeVector(this.x / L, this.y / L, this.z / L);
    },
    copyFrom(v) {
      this.x = v.x;
      this.y = v.y;
      this.z = v.z;
    },
    clone() {
      return makeVector(this.x, this.y, this.z);
    },
  };
}

// Only add GUI and methods if not already present (jest.setup may have partial BABYLON).
if (!global.BABYLON) global.BABYLON = {};
const B = global.BABYLON;

// Core vector & axis helpers (if missing)
if (!B.Vector3)
  B.Vector3 = function (x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  };
if (!B.Axis) B.Axis = { Z: makeVector(0, 0, 1), Y: makeVector(0, 1, 0) };

// Mesh & plane creation stubs
if (!B.MeshBuilder) B.MeshBuilder = {};
if (!B.MeshBuilder.CreatePlane) {
  B.MeshBuilder.CreatePlane = (name, opts, scene) => {
    return { name, opts, scene, position: makeVector(0, 0, 0), rotation: { y: 0 }, metadata: null };
  };
}

// GUI stubs sufficient for construction logic
if (!B.GUI) B.GUI = {};
if (!B.GUI.Control) {
  B.GUI.Control = { HORIZONTAL_ALIGNMENT_CENTER: 0 };
}
if (!B.GUI.AdvancedDynamicTexture) {
  B.GUI.AdvancedDynamicTexture = class {
    constructor() {
      this._rootContainer = { children: [] };
    }
    static CreateForMesh(_mesh, _w, _h) {
      return new B.GUI.AdvancedDynamicTexture();
    }
    addControl(c) {
      this._rootContainer.children.push(c);
    }
  };
}
if (!B.GUI.Button) {
  B.GUI.Button = class {
    static CreateSimpleButton(id, text) {
      const obs = () => ({ add: () => {} });
      return {
        id,
        text,
        width: '',
        height: '',
        color: '',
        thickness: 0,
        background: '',
        fontSize: 0,
        onPointerDownObservable: { add: () => {} },
        onPointerUpObservable: { add: () => {} },
        onPointerEnterObservable: { add: () => {} },
        onPointerOutObservable: { add: () => {} },
      };
    }
  };
}
if (!B.GUI.Rectangle) {
  B.GUI.Rectangle = class {
    constructor() {
      this.children = [];
      this.thickness = 0;
      this.cornerRadius = 0;
      this.background = '';
    }
    addControl(c) {
      this.children.push(c);
    }
  };
}
if (!B.GUI.StackPanel) {
  B.GUI.StackPanel = class {
    constructor() {
      this.isVertical = false;
      this.children = [];
    }
    addControl(c) {
      this.children.push(c);
    }
  };
}
if (!B.GUI.Image) {
  B.GUI.Image = class {
    constructor(id, source) {
      this.id = id;
      this.source = source;
      this.width = '';
      this.height = '';
      this.stretch = 0;
    }
    markAsDirty() {}
  };
  B.GUI.Image.STRETCH_UNIFORM = 0;
}
if (!B.GUI.TextBlock) {
  B.GUI.TextBlock = class {
    constructor(id, text) {
      this.id = id;
      this.text = text;
      this.color = '';
      this.fontSize = 0;
    }
  };
}

// Scene & camera stub (prefer existing Scene if defined)
const scene = new (B.Scene ||
  function () {
    this.onBeforeRenderObservable = {
      _l: [],
      add(fn) {
        this._l.push(fn);
        return fn;
      },
    };
  })();
scene.activeCamera = {
  fov: Math.PI / 2,
  position: makeVector(0, 0, 0),
  getDirection(axis) {
    if (axis === B.Axis.Z) return makeVector(0, 0, 1);
    if (axis === B.Axis.Y) return makeVector(0, 1, 0);
    return makeVector(0, 0, 0);
  },
};

// Provide simple observable add/remove if missing
if (!scene.onBeforeRenderObservable) {
  scene.onBeforeRenderObservable = {
    _l: [],
    add(fn) {
      this._l.push(fn);
      return fn;
    },
    remove(fn) {
      const i = this._l.indexOf(fn);
      if (i >= 0) this._l.splice(i, 1);
    },
  };
}

describe('XR HUD vertical positioning', () => {
  test('top energy plot lowered near center and bottom bar slightly raised', () => {
    const result = ensureWorldHUD({ scene });
    expect(result).toBeTruthy();
    const halfHeight = Math.tan((scene.activeCamera.fov || Math.PI / 2) / 2) * 1.1; // matches module logic

    // Simulate one frame so follow observers adjust (if they rely on onBeforeRender)
    for (const fn of scene.onBeforeRenderObservable._l) {
      fn();
    }

    const topY = result.planeTop.position.y; // should be small positive
    const bottomY = result.plane.position.y; // should be negative

    // Desired band selection (post-fix expectations):
    // Top: between 10% and 30% of halfHeight (legacy 60% is too high)
    const topMin = halfHeight * 0.1;
    const topMax = halfHeight * 0.3;
    expect(topY).toBeGreaterThanOrEqual(topMin);
    expect(topY).toBeLessThanOrEqual(topMax);

    // Bottom: magnitude between 35% and 55% (legacy 60% slightly farther; we want it a bit closer)
    const bottomAbs = Math.abs(bottomY);
    const bottomMin = halfHeight * 0.35;
    const bottomMax = halfHeight * 0.55;
    expect(bottomAbs).toBeGreaterThanOrEqual(bottomMin);
    expect(bottomAbs).toBeLessThanOrEqual(bottomMax);
  });
});
