import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.ts';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Minimal BABYLON scaffolding to satisfy moleculeView highlight path
if (!global.BABYLON) global.BABYLON = {};

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
  normalize() {
    const L = this.length() || 1;
    this.x /= L;
    this.y /= L;
    this.z /= L;
    return this;
  }
  normalizeToNew() {
    const L = this.length() || 1;
    return new Vector3(this.x / L, this.y / L, this.z / L);
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static Cross(a, b) {
    return new Vector3(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x,
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
  static FromEulerAngles(x, y, z) {
    const cx = Math.cos(x / 2);
    const sx = Math.sin(x / 2);
    const cy = Math.cos(y / 2);
    const sy = Math.sin(y / 2);
    const cz = Math.cos(z / 2);
    const sz = Math.sin(z / 2);
    return new Quaternion(
      sx * cy * cz - cx * sy * sz,
      cx * sy * cz + sx * cy * sz,
      cx * cy * sz - sx * sy * cz,
      cx * cy * cz + sx * sy * sz,
    );
  }
  static RotationAxis(axis, angle) {
    const half = angle / 2;
    const s = Math.sin(half);
    return new Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(half));
  }
}

class StandardMaterial {
  constructor() {
    this.diffuseColor = {};
    this.emissiveColor = {};
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
    position: Vector3.Zero(),
    scaling: new Vector3(1, 1, 1),
    rotationQuaternion: Quaternion.Identity(),
    thinInstanceSetBuffer() {},
    isPickable: true,
  };
}

const MeshBuilder = {
  CreateSphere: (name) => makeMaster(name),
  CreateCylinder: (name) => makeMaster(name),
  CreateLines: () => ({}),
};

Object.assign(global.BABYLON, {
  Vector3,
  Quaternion,
  StandardMaterial,
  Matrix,
  MeshBuilder,
  Axis: { X: new Vector3(1, 0, 0), Y: new Vector3(0, 1, 0), Z: new Vector3(0, 0, 1) },
});
global.BABYLON.Vector3.Up = () => new Vector3(0, 1, 0);

function makeScene() {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    getMeshByName() {
      return null;
    },
  };
}

describe('x-vr atom highlight alignment', () => {
  test('highlight remains centered on selected atom after rotation', () => {
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
    const selection = createSelectionService(state);
    selection.clickAtom(1);

    const highlight = view?._internals?.highlight?.atom;
    expect(highlight?.isVisible).toBe(true);
    expect(highlight.position.x).toBeCloseTo(1);

    const rot = Quaternion.FromEulerAngles(Math.PI / 4, Math.PI / 6, Math.PI / 8);
    for (const group of view._internals.atomGroups.values()) {
      group.master.rotationQuaternion = rot;
    }

    expect(highlight.position.x).toBeCloseTo(1, 4);
    expect(highlight.position.y).toBeCloseTo(0, 4);
    expect(highlight.position.z).toBeCloseTo(0, 4);
  });
});
