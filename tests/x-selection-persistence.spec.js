import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

// Reuse BABYLON stubs if not present (shared with other highlight tests).
if (!global.BABYLON) global.BABYLON = {};
const BAB = global.BABYLON;

BAB.Vector3 ||= class {
  constructor(x = 0, y = 0, z = 0) {
    this.x = x;
    this.y = y;
    this.z = z;
  }
  clone() {
    return new BAB.Vector3(this.x, this.y, this.z);
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
    return new BAB.Vector3(0, 0, 0);
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static Cross(a, b) {
    return new BAB.Vector3(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x
    );
  }
};

BAB.Vector3.Up = BAB.Vector3.Up || (() => new BAB.Vector3(0, 1, 0));

BAB.Quaternion ||= class {
  constructor(x = 0, y = 0, z = 0, w = 1) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
  }
  static Identity() {
    return new BAB.Quaternion();
  }
  static RotationAxis(axis, angle) {
    const half = angle / 2;
    const s = Math.sin(half);
    return new BAB.Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(half));
  }
};

BAB.Color3 ||= class {
  constructor(r = 0, g = 0, b = 0) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  clone() {
    return new BAB.Color3(this.r, this.g, this.b);
  }
  scale(f) {
    return new BAB.Color3(this.r * f, this.g * f, this.b * f);
  }
};

BAB.StandardMaterial ||= class {
  constructor() {
    this.diffuseColor = new BAB.Color3(1, 1, 1);
    this.emissiveColor = new BAB.Color3(0, 0, 0);
    this.specularColor = new BAB.Color3(0, 0, 0);
  }
};

BAB.Matrix ||= class {
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
    return new BAB.Matrix(arr);
  }
};

function makeMesh(name) {
  return {
    name,
    position: BAB.Vector3.Zero(),
    scaling: new BAB.Vector3(1, 1, 1),
    rotationQuaternion: BAB.Quaternion.Identity(),
    thinInstanceSetBuffer() {},
    thinInstanceEnablePicking: false,
    isPickable: false,
    parent: null,
    setEnabled() {},
  };
}

BAB.MeshBuilder ||= {
  CreateSphere: () => makeMesh('sphere'),
  CreateCylinder: () => makeMesh('cyl'),
  CreateLines: () => makeMesh('lines'),
};

BAB.Axis ||= { X: new BAB.Vector3(1, 0, 0), Y: new BAB.Vector3(0, 1, 0), Z: new BAB.Vector3(0, 0, 1) };

function makeScene() {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    getMeshByName() {
      return null;
    },
  };
}

function buildMolecule() {
  return {
    elements: [6, 6, 6],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1.5, y: 0, z: 0 },
      { x: 3.0, y: 0, z: 0 },
    ],
    bonds: [
      { i: 0, j: 1 },
      { i: 1, j: 2 },
    ],
  };
}

describe('x-selection-persistence', () => {
  test('bond selection persists through rotation operations', () => {
    const mol = buildMolecule();
    const state = createMoleculeState(mol);
    const selection = createSelectionService(state);
    const view = createMoleculeView(makeScene(), state);
    const manip = createManipulationService(state);

    selection.clickBond({ i: 0, j: 1, index: 0, key: '0-1' });
    expect(state.selection?.kind).toBe('bond');
    expect(view._internals.highlight.bond.isVisible).toBe(true);

    for (const angle of [5, 10, 15]) {
      const before = state.versions.positions;
      manip.rotateBond((angle * Math.PI) / 180);
      expect(state.versions.positions).toBeGreaterThanOrEqual(before);
      expect(state.selection?.kind).toBe('bond');
      expect(view._internals.highlight.bond.isVisible).toBe(true);
    }
  });

  test('atom selection remains after bondsChanged recompute', () => {
    const mol = buildMolecule();
    const state = createMoleculeState(mol);
    const selection = createSelectionService(state);
    const view = createMoleculeView(makeScene(), state);

    selection.clickAtom(1);
    expect(view._internals.highlight.atom.isVisible).toBe(true);
    state.markBondsChanged();
    expect(selection.get().kind).toBe('atom');
    expect(view._internals.highlight.atom.isVisible).toBe(true);
  });
});
