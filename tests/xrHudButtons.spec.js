/**
 * xrHudButtons.spec
 * Ensures PBC / Forces buttons use static labels and highlight (green) only when active,
 * and that bottom HUD bar height reduced (~30%).
 */

// Minimal BABYLON GUI stubs similar to existing positioning test (define before import)
if (!global.BABYLON) global.BABYLON = {};
/**
 * xrHudButtons.spec
 * Validates static labels + highlight background colors for PBC/Forces buttons
 * and reduced bottom bar height.
 */

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

// Reuse robust stubs from xrHudPosition test
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
    }
  };
}

import { ensureWorldHUD } from '../public/vr/xr-hud-bars.js';

describe('XR HUD buttons styling', () => {
  test('PBC & Forces static labels + highlight logic + reduced height', async () => {
    const scene = {
      activeCamera: {
        fov: Math.PI / 2,
        position: makeVector(0, 0, 0),
        getDirection(axis) {
          if (axis === B.Axis.Z) return makeVector(0, 0, 1);
          if (axis === B.Axis.Y) return makeVector(0, 1, 0);
          return makeVector(0, 0, 1);
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
    const viewer = {
      state: {
        showForces: false,
        showCell: false,
        cell: { enabled: false },
        toggleForceVectorsVisibility() {
          this.showForces = !this.showForces;
        },
        toggleCellVisibility() {
          this.showCell = !this.showCell;
          this.cell.enabled = this.showCell;
        },
        toggleGhostCells() {},
      },
    };
    window._viewer = viewer;
    const res = ensureWorldHUD({ scene, getViewer: () => viewer });
    expect(res).toBeTruthy();
    expect(res.plane.opts.height).toBeGreaterThan(0.25);
    expect(res.plane.opts.height).toBeLessThan(0.28);
    const fb = window.__XR_HUD_FALLBACK;
    expect(fb).toBeTruthy();
    function findBtn(label) {
      const rootChildren = fb.tex?._rootContainer?.children || [];
      const q = [...rootChildren];
      while (q.length) {
        const c = q.shift();
        if (!c) continue;
        if (c.text === label || c?.textBlock?.text === label) return c;
        if (Array.isArray(c.children)) q.push(...c.children);
      }
      return null;
    }
    const btnForces = findBtn('Forces');
    const beforeForces = btnForces?.background;
    // Wait for async controlsModel import hookup (poll a few ticks)
    let cm = null;
    for (let i = 0; i < 10 && !cm; i++) {
      cm = res.getControlsModel && res.getControlsModel();
      if (!cm) await new Promise((r) => setTimeout(r, 0));
    }
    expect(cm).toBeTruthy();
    cm.toggleForces();
    // Allow onStateChange -> syncButtons to run (microtask)
    await new Promise((r) => setTimeout(r, 0));
    // Extra frames (though not strictly needed now)
    for (let i = 0; i < 2; i++) scene.onBeforeRenderObservable._l.forEach((fn) => fn());
    const afterForces = btnForces?.background;
    expect(afterForces).not.toBe(beforeForces);
    expect(afterForces).toBe('rgba(40,140,100,0.9)');
    const btnPbc = findBtn('PBC');
    const beforePbc = btnPbc?.background;
    cm.togglePBC();
    await new Promise((r) => setTimeout(r, 0));
    for (let i = 0; i < 2; i++) scene.onBeforeRenderObservable._l.forEach((fn) => fn());
    const afterPbc = btnPbc?.background;
    expect(afterPbc).not.toBe(beforePbc);
    expect(afterPbc).toBe('rgba(40,140,100,0.9)');
  });
});
// (End of file)
