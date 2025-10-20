/**
 * xrHudEnergyPlotPicking.spec
 * Validates that the top energy plot HUD panel:
 *  - Is non-pickable by default (so it doesn't block atom/bond selection rays)
 *  - Is positioned at ~2x distance compared to bottom bar when using defaults
 *  - Is scaled up accordingly (width & height ~2x)
 * Provides an override test ensuring options can restore pickable behavior and base distance.
 */

import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

// Provide window/document & canvas stub (energy plot requires toDataURL)
if (!global.window) global.window = { location: { search: '' } };
if (!global.document) {
  global.document = {
    createElement: (tag) =>
      tag === 'canvas'
        ? { width: 0, height: 0, getContext: () => ({}), toDataURL: () => 'data:image/png;base64,' }
        : { style: {} },
  };
}

// BABYLON stubs (similar to other HUD tests)
if (!global.BABYLON) global.BABYLON = {};
const B = global.BABYLON;
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
if (!B.Vector3)
  B.Vector3 = function (x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  };
if (!B.Axis) B.Axis = { Z: makeVector(0, 0, 1), Y: makeVector(0, 1, 0) };
if (!B.MeshBuilder) B.MeshBuilder = {};
if (!B.MeshBuilder.CreatePlane) {
  B.MeshBuilder.CreatePlane = (name, opts, scene) => ({
    name,
    opts,
    scene,
    position: makeVector(0, 0, 0),
    rotation: { y: 0 },
    metadata: null,
    isPickable: true,
  });
}
if (!B.GUI) B.GUI = {};
if (!B.GUI.Control)
  B.GUI.Control = { HORIZONTAL_ALIGNMENT_CENTER: 0, HORIZONTAL_ALIGNMENT_LEFT: 1 };
if (!B.GUI.AdvancedDynamicTexture) {
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
}
if (!B.GUI.Button) {
  B.GUI.Button = class {
    static CreateSimpleButton(id, text) {
      return {
        id,
        text,
        textBlock: { text },
        width: '',
        height: '',
        background: '',
        thickness: 0,
        fontSize: 0,
        children: [],
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
      this.background = '';
      this.cornerRadius = 0;
      this.thickness = 0;
    }
    addControl(c) {
      this.children.push(c);
    }
  };
}
if (!B.GUI.StackPanel) {
  B.GUI.StackPanel = class {
    constructor() {
      this.children = [];
      this.isVertical = false;
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
      this.textHorizontalAlignment = 0;
    }
  };
}

function makeScene() {
  const scene = {
    onBeforeRenderObservable: {
      _l: [],
      add(fn) {
        this._l.push(fn);
        return fn;
      },
    },
    activeCamera: {
      fov: Math.PI / 2,
      position: makeVector(0, 0, 0),
      getDirection(axis) {
        if (axis === B.Axis.Z) return makeVector(0, 0, 1);
        if (axis === B.Axis.Y) return makeVector(0, 1, 0);
        return makeVector(0, 0, 0);
      },
    },
  };
  return scene;
}

function runFrames(scene, n = 1) {
  for (let i = 0; i < n; i++) scene.onBeforeRenderObservable._l.forEach((fn) => fn());
}

describe('XR HUD energy plot picking & distance', () => {
  test('default: top energy plane non-pickable, at 2x distance & scaled', () => {
    const scene = makeScene();
    const res = ensureWorldHUD({ scene });
    expect(res).toBeTruthy();
    runFrames(scene, 2); // allow follow observers to update positions
    const bottomZ = res.plane.position.z;
    const topZ = res.planeTop.position.z;
    // Base distance is 1.1 -> expect bottom ~1.1 and top ~2.2 (allow tolerance)
    expect(bottomZ).toBeGreaterThan(1.0);
    expect(bottomZ).toBeLessThan(1.2);
    expect(topZ).toBeGreaterThan(2.0);
    expect(topZ).toBeLessThan(2.4);
    // Non-pickable default
    expect(res.planeTop.isPickable).toBe(false);
    // Size ratio reflects default internal plane width difference times scale multiplier.
    // Bottom width = 2.0, top base width = 1.35*4 = 5.4, then scaled by topEnergyScaleMult (default 2) => 10.8 => ratio 5.4.
    const widthRatio = res.planeTop.opts.width / res.plane.opts.width;
    expect(widthRatio).toBeGreaterThan(5.3);
    expect(widthRatio).toBeLessThan(5.5);
  });

  test('overrides: pickable true & base distance restored', () => {
    const scene = makeScene();
    // Clear previous cached HUD if any
    if (global.window && window.__XR_HUD_FALLBACK) {
      delete window.__XR_HUD_FALLBACK;
      delete window.__XR_HUD_FALLBACK_TOP;
    }
    const res = ensureWorldHUD({
      scene,
      topEnergyDistanceMult: 1,
      topEnergyScaleMult: 1,
      topEnergyIsPickable: true,
    });
    expect(res).toBeTruthy();
    runFrames(scene, 2);
    const bottomZ = res.plane.position.z;
    const topZ = res.planeTop.position.z;
    // Distances roughly equal now (since multiplier = 1)
    expect(Math.abs(topZ - bottomZ)).toBeLessThan(0.1);
    expect(res.planeTop.isPickable).toBe(true);
    // Width ratio falls back to original (without extra scaleMult) -> base top width 5.4 / bottom 2.0 = 2.7
    const widthRatio = res.planeTop.opts.width / res.plane.opts.width;
    expect(widthRatio).toBeGreaterThan(2.6);
    expect(widthRatio).toBeLessThan(2.8);
  });
});
