/** @jest-environment jsdom */

// Verify that force recomputes refresh the thin-instance buffers backing force vectors.

if (!global.BABYLON) global.BABYLON = {};
const BAB = global.BABYLON;

class Vector3 {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  clone() {
    return new Vector3(this.x, this.y, this.z);
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
    return this.clone().normalize();
  }
  static Zero() {
    return new Vector3(0, 0, 0);
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static Cross(a, b) {
    return new Vector3(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x
    );
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
    const half = angle / 2;
    const s = Math.sin(half);
    return new Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(half));
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

class Matrix {
  constructor(values) {
    this.m = values || new Float32Array(16);
  }
  static Compose(scale, rot, pos) {
    const arr = new Float32Array(16);
    arr[0] = scale.x;
    arr[5] = scale.y;
    arr[10] = scale.z;
    arr[12] = pos.x;
    arr[13] = pos.y;
    arr[14] = pos.z;
    return new Matrix(arr);
  }
}

class StandardMaterial {
  constructor(name) {
    this.name = name;
    this.diffuseColor = new Color3(1, 1, 1);
    this.emissiveColor = new Color3(0, 0, 0);
    this.specularColor = new Color3(0, 0, 0);
    this.disableLighting = false;
  }
}

class FakeMesh {
  constructor(name) {
    this.name = name;
    this._buffers = {};
    this.isPickable = false;
    this.thinInstanceEnablePicking = false;
    this.material = null;
    this.parent = null;
    this.isVisible = false;
  }
  thinInstanceSetBuffer(kind, data) {
    this._buffers[kind] = data;
  }
  thinInstanceRefreshBoundingInfo() {}
  setEnabled(on) {
    this.isVisible = !!on;
  }
}

const MeshBuilder = {
  CreateCylinder: (name) => new FakeMesh(name),
  CreateSphere: (name) => new FakeMesh(name),
  CreateLines: () => new FakeMesh('lines'),
};

BAB.Vector3 ||= Vector3;
BAB.Quaternion ||= Quaternion;
BAB.Color3 ||= Color3;
BAB.StandardMaterial ||= StandardMaterial;
BAB.Matrix ||= Matrix;
BAB.MeshBuilder ||= MeshBuilder;
BAB.Axis ||= { X: new Vector3(1, 0, 0), Y: new Vector3(0, 1, 0), Z: new Vector3(0, 0, 1) };
BAB.Vector3.Up = () => new Vector3(0, 1, 0);

jest.mock('../public/render/scene.js', () => ({
  createScene: async () => ({
    engine: { runRenderLoop: () => {}, stopRenderLoop: () => {} },
    scene: {
      meshes: [],
      render: () => {},
      dispose: () => {},
      onPointerObservable: { add() {} },
      onBeforeRenderObservable: { add() {} },
      getMeshByName: () => null,
    },
    camera: { attachControl: () => {}, detachControl: () => {} },
  }),
}));

const firstForces = [
  [0.0, 0.1, 0.0],
  [0.05, 0.0, 0.0],
  [-0.02, 0.03, 0.0],
];
const secondForces = [
  [0.5, 0.0, 0.0],
  [0.1, 0.1, 0.0],
  [-0.3, 0.2, 0.0],
];

describe('x-force-update-perturbation', () => {
  let origFetch;
  let origWebSocket;

  beforeAll(() => {
    origFetch = global.fetch;
    origWebSocket = global.WebSocket;
    let calls = 0;
    global.fetch = jest.fn(async () => {
      calls += 1;
      const payload =
        calls <= 1
          ? { results: { energy: -10.0, forces: firstForces } }
          : { results: { energy: -12.0, forces: secondForces } };
      return { ok: true, status: 200, json: async () => payload };
    });
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
      onmessage() {}
      onerror() {}
    }
    global.WebSocket = FakeWS;
  });

  afterAll(() => {
    global.fetch = origFetch;
    global.WebSocket = origWebSocket;
  });

  test('force vector instance buffer changes after recompute', async () => {
    window.__MLIPVIEW_TEST_MODE = true;
    const canvas = document.createElement('canvas');
    canvas.id = 'viewer';
    canvas.width = 800;
    canvas.height = 600;
    canvas.getContext = () => ({
      clearRect() {},
      beginPath() {},
      moveTo() {},
      lineTo() {},
      stroke() {},
      arc() {},
      fill() {},
    });
    document.body.appendChild(canvas);

    const { initNewViewer } = await import('../public/index.js');
    const viewer = await initNewViewer(canvas, {
      elements: ['O', 'H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 0.96, y: 0, z: 0 },
        { x: -0.24, y: 0.93, z: 0 },
      ],
      bonds: [],
    });

    viewer.setForceVectorsEnabled(true);
    viewer.state.forces = firstForces.map((f) => f.slice());
    viewer.state.bus.emit('forcesChanged');
    const group = viewer.view._internals.forceGroups.get('force');
    expect(group).toBeTruthy();
    const captureMatrices = () =>
      group.mats.map((mat) => (mat?.m ? Array.from(mat.m) : []));

    let initialMats = captureMatrices();
    for (let i = 0; (!initialMats.length || !initialMats[0].length) && i < 20; i++) {
      viewer.view.rebuildForces?.();
      // eslint-disable-next-line no-await-in-loop
      await new Promise((r) => setTimeout(r, 5));
      initialMats = captureMatrices();
    }
    expect(initialMats.length).toBeGreaterThan(0);
    expect(initialMats[0].length).toBe(16);

    viewer.state.positions[0].x += 5;
    viewer.state.positions[0].y += 4;
    viewer.state.positions[0].z -= 3;
    viewer.state.markPositionsChanged();

    viewer.state.forces = secondForces.map((f) => f.slice());
    viewer.state.bus.emit('forcesChanged');
    let afterMats = captureMatrices();
    for (let i = 0; (!afterMats.length || !afterMats[0].length) && i < 20; i++) {
      viewer.view.rebuildForces?.();
      // eslint-disable-next-line no-await-in-loop
      await new Promise((r) => setTimeout(r, 5));
      afterMats = captureMatrices();
    }
    expect(afterMats.length).toBeGreaterThan(0);
    let changed = false;
    for (let idx = 0; idx < Math.min(initialMats.length, afterMats.length); idx++) {
      const before = initialMats[idx];
      const after = afterMats[idx];
      if (before.length !== after.length) {
        changed = true;
        break;
      }
      for (let j = 0; j < before.length; j++) {
        if (before[j] !== after[j]) {
          changed = true;
          break;
        }
      }
      if (changed) break;
    }
    expect(changed).toBe(true);
    expect(Array.isArray(viewer.state.forces)).toBe(true);
    expect(viewer.state.forces.length).toBe(3);
  });
});
