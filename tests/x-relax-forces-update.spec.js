/** @jest-environment jsdom */

import { getWS } from '../public/fairchem_ws_client.js';

if (!global.BABYLON) global.BABYLON = {};
const BABYLON = global.BABYLON;

class Vector3 {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  length() {
    return Math.hypot(this.x, this.y, this.z);
  }
  normalize() {
    const L = this.length() || 1;
    this.x /= L;
    this.y /= L;
    this.z /= L;
    return this;
  }
  normalizeToNew() {
    return new Vector3(this.x, this.y, this.z).normalize();
  }
  static Cross(a, b) {
    return new Vector3(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x
    );
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
}

class Quaternion {
  constructor(x = 0, y = 0, z = 0, w = 1) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
  }
  static Identity() {
    return new Quaternion();
  }
  static RotationAxis(axis, angle) {
    const s = Math.sin(angle / 2);
    return new Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(angle / 2));
  }
}

class Matrix {
  constructor() {
    this.m = new Float32Array(16);
  }
  static Compose(scale, _rot, pos) {
    const mat = new Matrix();
    mat.m[12] = pos.x;
    mat.m[13] = pos.y;
    mat.m[14] = pos.z;
    return mat;
  }
}

class Color3 {
  constructor(r = 0, g = 0, b = 0) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  clone() {
    return new Color3(this.r, this.g, this.b);
  }
  scale(f) {
    return new Color3(this.r * f, this.g * f, this.b * f);
  }
}

class StandardMaterial {
  constructor() {
    this.diffuseColor = new Color3(1, 1, 1);
    this.emissiveColor = new Color3();
    this.specularColor = new Color3();
    this.alpha = 1;
  }
}

class Mesh {
  constructor() {
    this._buffers = {};
    this.isVisible = true;
  }
  thinInstanceSetBuffer(kind, arr) {
    this._buffers[kind] = arr;
  }
  setEnabled(on) {
    this.isVisible = on;
  }
}

const MeshBuilder = {
  CreateSphere: () => new Mesh(),
  CreateCylinder: () => new Mesh(),
  CreateLines: () => new Mesh(),
};

class TransformNode {
  constructor() {
    this.position = new Vector3();
    this.scaling = new Vector3(1, 1, 1);
    this.parent = null;
  }
}

Object.assign(BABYLON, {
  Vector3,
  Quaternion,
  Matrix,
  Color3,
  StandardMaterial,
  MeshBuilder,
  TransformNode,
  Material: { MATERIAL_ALPHABLEND: 2 },
  PointerEventTypes: { POINTERDOWN: 1 },
});

window.__MLIPVIEW_TEST_MODE = true;
window.FORCE_MIN = 0.01;

jest.mock('../public/render/scene.js', () => {
  class mockVec3 {
    constructor(x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
  }
  class mockEngine {
    runRenderLoop(cb) {
      for (let i = 0; i < 2; i++) cb();
    }
    stopRenderLoop() {}
  }
  class mockScene {
    constructor() {
      this.onBeforeRenderObservable = { add() {} };
      this.onPointerObservable = { add() {} };
      this.pointerX = 0;
      this.pointerY = 0;
    }
    render() {}
    getEngine() {
      return { getRenderingCanvas: () => ({ addEventListener() {} }) };
    }
    createPickingRay() {
      return { origin: new mockVec3(), direction: new mockVec3(0, 1, 0) };
    }
    pick() {
      return { hit: false };
    }
  }
  return {
    createScene: async () => ({
      engine: new mockEngine(),
      scene: new mockScene(),
      camera: {
        detachControl() {},
        attachControl() {},
      },
    }),
  };
});

const baseForces = [
  [0.1, 0.0, 0.0],
  [0.0, 0.1, 0.0],
  [0.0, 0.0, 0.1],
];

function scaledForces(step) {
  return baseForces.map((f) => [
    f[0] * (1 + 0.2 * step),
    f[1] * (1 + 0.15 * step),
    f[2] * (1 + 0.1 * step),
  ]);
}

function forceInstanceCount(api) {
  const group = api.view._internals.forceGroups.get('force');
  if (!group) return 0;
  const buf = group.master._buffers['matrix'];
  return buf ? buf.length / 16 : 0;
}

describe('x-relax-forces-update', () => {
  test('multiple relax steps update force vector buffers', async () => {
    const origWS = global.WebSocket;
    class FakeWS {
      constructor() {
        this.readyState = 0;
        setTimeout(() => {
          this.readyState = 1;
          this.onopen?.();
        }, 0);
      }
      send() {}
      close() {}
      onopen() {}
      onerror() {}
    }
    global.WebSocket = FakeWS;
    const ws = getWS();
    ws.setTestHook(() => {});

    const { initNewViewer } = await import('../public/index.js');
    const canvas = document.createElement('canvas');
    canvas.id = 'viewer';
    canvas.addEventListener = () => {};
    document.body.appendChild(canvas);

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
      fillRect() {},
    });
    energyWrapper.appendChild(energyCanvas);
    energyWrapper.appendChild(document.createElement('div')).id = 'energyLabel';

    const api = await initNewViewer(canvas, {
      elements: ['O', 'H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.95, y: 0, z: 0 },
        { x: -0.24, y: 0.93, z: 0 },
      ],
      bonds: [],
    });

    await new Promise((r) => setTimeout(r, 0));
    ws.injectTestResult({
      positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
      forces: scaledForces(0),
      energy: -5,
    });
    api.state.bus.emit('forcesChanged');
    for (let i = 0; i < 20 && forceInstanceCount(api) === 0; i++) {
      api.view.rebuildForces();
      // eslint-disable-next-line no-await-in-loop
      await new Promise((r) => setTimeout(r, 5));
    }

    const initialMatrix = api.view._internals.forceGroups
      .get('force')
      .master._buffers['matrix'].slice();

    for (let step = 1; step <= 10; step++) {
      const relaxPromise = api.relaxStep();
      await new Promise((r) => setTimeout(r, 0));
      ws.injectTestResult({
        positions: api.state.positions.map((p) => [p.x, p.y, p.z]),
        forces: scaledForces(step),
        energy: -5 - 0.1 * step,
      });
      await relaxPromise;
      api.state.bus.emit('forcesChanged');
      api.view.rebuildForces();
      // eslint-disable-next-line no-await-in-loop
      await new Promise((r) => setTimeout(r, 5));
    }

    const afterMatrix = api.view._internals.forceGroups.get('force').master._buffers['matrix'];
    let changed = false;
    for (let i = 0; i < afterMatrix.length && i < initialMatrix.length; i++) {
      if (afterMatrix[i] !== initialMatrix[i]) {
        changed = true;
        break;
      }
    }
    expect(forceInstanceCount(api)).toBeGreaterThan(0);
    expect(changed).toBe(true);

    ws.setTestHook(null);
    global.WebSocket = origWS;
  });
});
