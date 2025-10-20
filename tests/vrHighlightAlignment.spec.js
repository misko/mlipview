// VR highlight alignment regression test
// Ensures that after applying a rotation (simulating VR controller induced rotation)
// the cyan atom highlight sphere remains centered on the selected atom.

import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Provide minimal BABYLON implementation sufficient for rotation + parenting behavior
global.BABYLON = global.BABYLON || {};

class Vector3 {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  clone() {
    return new Vector3(this.x, this.y, this.z);
  }
  static Zero() {
    return new Vector3(0, 0, 0);
  }
  length() {
    return Math.hypot(this.x, this.y, this.z);
  }
  normalizeToNew() {
    const L = this.length() || 1;
    return new Vector3(this.x / L, this.y / L, this.z / L);
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static Cross(a, b) {
    return new Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
  }
  normalize() {
    const L = this.length() || 1;
    this.x /= L;
    this.y /= L;
    this.z /= L;
    return this;
  }
  add(v) {
    this.x += v.x;
    this.y += v.y;
    this.z += v.z;
    return this;
  }
  scale(f) {
    this.x *= f;
    this.y *= f;
    this.z *= f;
    return this;
  }
  subtract(v) {
    return new Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
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
  static FromEulerAngles(x, y, z) {
    // Simplified ZYX intrinsic
    const cx = Math.cos(x / 2),
      sx = Math.sin(x / 2);
    const cy = Math.cos(y / 2),
      sy = Math.sin(y / 2);
    const cz = Math.cos(z / 2),
      sz = Math.sin(z / 2);
    return new Quaternion(
      sx * cy * cz - cx * sy * sz,
      cx * sy * cz + sx * cy * sz,
      cx * cy * sz - sx * sy * cz,
      cx * cy * cz + sx * sy * sz
    );
  }
  clone() {
    return new Quaternion(this.x, this.y, this.z, this.w);
  }
  static RotationAxis(axis, angle) {
    const half = angle / 2;
    const s = Math.sin(half);
    const c = Math.cos(half);
    const ax = axis.x,
      ay = axis.y,
      az = axis.z;
    return new Quaternion(ax * s, ay * s, az * s, c);
  }
}
class StandardMaterial {
  constructor(name, scene) {
    this.name = name;
    this.diffuseColor = {};
    this.emissiveColor = {};
    this.alpha = 1;
  }
}
class Matrix {
  constructor() {
    this.m = new Float32Array(16);
  }
  static Compose() {
    return new Matrix();
  }
}
function makeMaster(name) {
  return {
    name,
    thinInstanceSetBuffer() {},
    position: Vector3.Zero(),
    scaling: new Vector3(1, 1, 1),
    rotationQuaternion: Quaternion.Identity(),
    isPickable: true,
  };
}
const MeshBuilder = {
  CreateSphere: (name) => makeMaster(name),
  CreateCylinder: (name) => makeMaster(name),
  CreateLines: () => ({}),
};

Object.assign(BABYLON, {
  Vector3,
  Quaternion,
  StandardMaterial,
  Matrix,
  MeshBuilder,
  Axis: { X: new Vector3(1, 0, 0), Y: new Vector3(0, 1, 0), Z: new Vector3(0, 0, 1) },
});

function makeScene() {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    getMeshByName() {
      return null;
    },
  };
}

describe('vr atom highlight alignment', () => {
  test('highlight sphere tracks atom after rotation', () => {
    const state = createMoleculeState({
      elements: ['C', 'H', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
        { x: 0, y: 1, z: 0 },
      ],
      bonds: [
        { i: 0, j: 1 },
        { i: 0, j: 2 },
      ],
    });
    const scene = makeScene();
    const view = createMoleculeView(scene, state);
    const sel = createSelectionService(state);
    sel.clickAtom(1); // select atom at (1,0,0)
    const hi = view._internals.highlight.atom;
    expect(hi.isVisible).toBe(true);
    // Record initial relative position before rotation
    const before = hi.position.clone();
    expect(before.x).toBeCloseTo(1, 5); // should be at the atom's local position
    // Simulate VR rotation by applying quaternion to all master atom meshes
    const rot = Quaternion.FromEulerAngles(Math.PI / 4, Math.PI / 6, 0); // arbitrary rotation
    for (const g of view._internals.atomGroups.values()) {
      g.master.rotationQuaternion = rot.clone();
    }
    // After rotation the highlight remains parented; its local position (1,0,0) unchanged
    const after = hi.position;
    expect(after.x).toBeCloseTo(1, 5);
    expect(after.y).toBeCloseTo(0, 5);
    expect(after.z).toBeCloseTo(0, 5);
  });
});
