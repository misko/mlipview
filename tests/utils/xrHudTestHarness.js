/**
 * Shared Babylon GUI stubs + helpers for xr-hud-bars unit tests.
 */

function vec(x = 0, y = 0, z = 0) {
  return {
    x,
    y,
    z,
    add(o) {
      return vec(this.x + o.x, this.y + o.y, this.z + o.z);
    },
    subtract(o) {
      return vec(this.x - o.x, this.y - o.y, this.z - o.z);
    },
    scale(s) {
      return vec(this.x * s, this.y * s, this.z * s);
    },
    normalize() {
      const L = Math.hypot(this.x, this.y, this.z) || 1;
      return vec(this.x / L, this.y / L, this.z / L);
    },
    copyFrom(o) {
      this.x = o.x;
      this.y = o.y;
      this.z = o.z;
      return this;
    },
    clone() {
      return vec(this.x, this.y, this.z);
    },
  };
}

function ensureDocumentCanvas() {
  if (global.document) return;
  global.document = {
    createElement(tag) {
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

function ensureWindow() {
  global.window = global.window || {};
  if (!window.location) window.location = { search: '' };
}

function ensureBabylon() {
  const B = (global.BABYLON = global.BABYLON || {});

  if (!B.Vector3)
    B.Vector3 = function (x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    };
  if (!B.Axis) B.Axis = { Z: vec(0, 0, 1), Y: vec(0, 1, 0) };
  if (!B.Mesh) B.Mesh = { BILLBOARDMODE_NONE: 0 };
  B.MeshBuilder = B.MeshBuilder || {};
  if (!B.MeshBuilder.CreatePlane) {
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
  }

  B.GUI = B.GUI || {};
  B.GUI.Control = B.GUI.Control || {
    HORIZONTAL_ALIGNMENT_CENTER: 0,
    HORIZONTAL_ALIGNMENT_LEFT: 1,
  };

  if (!B.GUI.AdvancedDynamicTexture) {
    B.GUI.AdvancedDynamicTexture = class AdvancedDynamicTexture {
      constructor(mesh, width, height) {
        this._rootContainer = { children: [] };
        this._linkedMesh = mesh || null;
        this.width = width || 1024;
        this.height = height || 512;
      }
      static CreateForMesh(mesh, width, height) {
        return new B.GUI.AdvancedDynamicTexture(mesh, width, height);
      }
      addControl(ctrl) {
        this._rootContainer.children.push(ctrl);
      }
    };
  }

  if (!B.GUI.StackPanel) {
    B.GUI.StackPanel = class StackPanel {
      constructor() {
        this.children = [];
        this.isVertical = true;
        this.horizontalAlignment = 0;
      }
      addControl(ctrl) {
        this.children.push(ctrl);
      }
    };
  }

  if (!B.GUI.Rectangle) {
    B.GUI.Rectangle = class Rectangle {
      constructor() {
        this.children = [];
        this.background = '';
        this.thickness = 0;
        this.cornerRadius = 0;
      }
      addControl(ctrl) {
        this.children.push(ctrl);
      }
    };
  }

  if (!B.GUI.TextBlock) {
    B.GUI.TextBlock = class TextBlock {
      constructor(id, text) {
        this.id = id;
        this.text = text || '';
        this.color = '#fff';
        this.fontSize = 16;
        this.textHorizontalAlignment = 0;
      }
    };
  }

  if (!B.GUI.Image) {
    B.GUI.Image = class Image {
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

  if (!B.GUI.Button) {
    B.GUI.Button = {
      CreateSimpleButton(id, text) {
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
          onPointerDownObservable: { add() {} },
          onPointerUpObservable: { add() {} },
          onPointerEnterObservable: { add() {} },
          onPointerOutObservable: { add() {} },
        };
      },
    };
  }
}

export function setupXRHudTestEnv() {
  ensureWindow();
  ensureDocumentCanvas();
  ensureBabylon();
  return global.BABYLON;
}

export function resetXRHudSingleton() {
  if (typeof window === 'undefined') return;
  delete window.__XR_HUD_FALLBACK;
  delete window.__XR_HUD_FALLBACK_TOP;
}

export function makeHudScene() {
  setupXRHudTestEnv();
  const observers = [];
  const scene = {
    activeCamera: {
      fov: Math.PI / 2,
      position: vec(),
      getDirection(axis) {
        if (axis === global.BABYLON.Axis.Z) return vec(0, 0, 1);
        if (axis === global.BABYLON.Axis.Y) return vec(0, 1, 0);
        return vec();
      },
    },
    onBeforeRenderObservable: {
      _l: observers,
      add(fn) {
        observers.push(fn);
        return fn;
      },
    },
    getEngine() {
      return {
        getCaps() {
          return { maxTextureSize: 2048 };
        },
      };
    },
  };
  return scene;
}

export function runHudFrames(scene, count = 1) {
  for (let i = 0; i < count; i += 1) {
    if (scene?.onBeforeRenderObservable?._l) {
      for (const fn of scene.onBeforeRenderObservable._l) {
        try {
          fn();
        } catch {}
      }
    }
  }
}

export function findGuiNode(root, predicate) {
  if (!root) return null;
  const q = [root];
  while (q.length) {
    const node = q.shift();
    if (!node) continue;
    if (predicate(node)) return node;
    if (Array.isArray(node.children)) q.push(...node.children);
  }
  return null;
}

export default {
  setupXRHudTestEnv,
  resetXRHudSingleton,
  makeHudScene,
  runHudFrames,
  findGuiNode,
};

