/**
 * xrHudEnergyDepth.spec
 * Verifies that increasing topEnergyDistanceMult pushes top panel forward (larger z) while Y band remains stable.
 */
import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';
if (!global.window) global.window = { location: { search: '' } };
if (!global.document)
  global.document = {
    createElement: (tag) =>
      tag === 'canvas' ? { width: 0, height: 0, getContext: () => ({}), toDataURL: () => '' } : {},
  };
if (!global.BABYLON) global.BABYLON = {};
const B = global.BABYLON;
function v(x = 0, y = 0, z = 0) {
  return {
    x,
    y,
    z,
    add(o) {
      return v(x + o.x, y + o.y, z + o.z);
    },
    subtract(o) {
      return v(x - o.x, y - o.y, z - o.z);
    },
    scale(s) {
      return v(x * s, y * s, z * s);
    },
    normalize() {
      const L = Math.hypot(x, y, z) || 1;
      return v(x / L, y / L, z / L);
    },
    copyFrom(o) {
      this.x = o.x;
      this.y = o.y;
      this.z = o.z;
    },
  };
}
if (!B.Axis) B.Axis = { Z: v(0, 0, 1), Y: v(0, 1, 0) };
if (!B.MeshBuilder) B.MeshBuilder = {};
if (!B.MeshBuilder.CreatePlane)
  B.MeshBuilder.CreatePlane = (n, o, s) => ({
    name: n,
    opts: o,
    scene: s,
    position: v(),
    rotation: { y: 0 },
    isPickable: true,
  });
if (!B.GUI) B.GUI = {};
if (!B.GUI.Control)
  B.GUI.Control = { HORIZONTAL_ALIGNMENT_CENTER: 0, HORIZONTAL_ALIGNMENT_LEFT: 1 };
if (!B.GUI.AdvancedDynamicTexture)
  B.GUI.AdvancedDynamicTexture = class {
    constructor() {
      this._rootContainer = { children: [] };
    }
    static CreateForMesh() {
      return new B.GUI.AdvancedDynamicTexture();
    }
    addControl(c) {
      this._rootContainer.children.push(c);
    }
  };
if (!B.GUI.StackPanel)
  B.GUI.StackPanel = class {
    constructor() {
      this.children = [];
    }
    addControl(c) {
      this.children.push(c);
    }
  };
if (!B.GUI.Rectangle)
  B.GUI.Rectangle = class {
    constructor() {
      this.children = [];
    }
    addControl(c) {
      this.children.push(c);
    }
  };
if (!B.GUI.Button)
  B.GUI.Button = class {
    static CreateSimpleButton(id, text) {
      return {
        id,
        text,
        background: '',
        onPointerDownObservable: { add: () => {} },
        onPointerUpObservable: { add: () => {} },
        onPointerEnterObservable: { add: () => {} },
        onPointerOutObservable: { add: () => {} },
      };
    }
  };
if (!B.GUI.Image) {
  B.GUI.Image = class {
    constructor(id, src) {
      this.id = id;
      this.source = src;
      this.width = '';
      this.height = '';
      this.stretch = 0;
    }
    markAsDirty() {}
  };
  B.GUI.Image.STRETCH_UNIFORM = 0;
}
if (!B.GUI.TextBlock)
  B.GUI.TextBlock = class {
    constructor(id, text) {
      this.id = id;
      this.text = text;
      this.color = '';
    }
  };

function makeScene() {
  return {
    activeCamera: {
      fov: Math.PI / 2,
      position: v(),
      getDirection(ax) {
        if (ax === B.Axis.Z) return v(0, 0, 1);
        if (ax === B.Axis.Y) return v(0, 1, 0);
        return v();
      },
    },
    onBeforeRenderObservable: {
      _l: [],
      add(fn) {
        this._l.push(fn);
        return fn;
      },
    },
  };
}

describe('XR HUD energy depth scaling', () => {
  test('topEnergyDistanceMult increases z without altering Y band drastically', () => {
    const scene1 = makeScene();
    if (window.__XR_HUD_FALLBACK) {
      delete window.__XR_HUD_FALLBACK;
      delete window.__XR_HUD_FALLBACK_TOP;
    }
    const res1 = ensureWorldHUD({ scene: scene1, topEnergyDistanceMult: 1 });
    scene1.onBeforeRenderObservable._l.forEach((fn) => fn());
    if (window.__XR_HUD_FALLBACK) {
      delete window.__XR_HUD_FALLBACK;
      delete window.__XR_HUD_FALLBACK_TOP;
    }
    const scene2 = makeScene();
    const res2 = ensureWorldHUD({ scene: scene2, topEnergyDistanceMult: 3 });
    scene2.onBeforeRenderObservable._l.forEach((fn) => fn());
    expect(res2.planeTop.position.z).toBeGreaterThan(res1.planeTop.position.z * 2.5); // 3x vs 1x (allow some tolerance)
    // Y should still be within reasonable factor (< 1.5x) of original since vertical uses base dist
    expect(Math.abs(res2.planeTop.position.y)).toBeLessThan(
      Math.abs(res1.planeTop.position.y) * 1.2 + 1e-6
    );
  });
});
