import fs from 'fs';
import path from 'path';
import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

global.BABYLON = global.BABYLON || {};
const BABYLON = global.BABYLON;

BABYLON.Vector3 ||= class {
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
    return new BABYLON.Vector3(this.x, this.y, this.z).normalize();
  }
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static Cross(a, b) {
    return new BABYLON.Vector3(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x
    );
  }
};

BABYLON.Quaternion ||= { Identity: () => ({}), RotationAxis: () => ({}) };
BABYLON.Matrix ||= class {
  constructor() {
    this.m = new Float32Array(16);
  }
  static Compose() {
    return new BABYLON.Matrix();
  }
};
BABYLON.Color3 ||= class {
  constructor(r = 1, g = 1, b = 1) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  clone() {
    return new BABYLON.Color3(this.r, this.g, this.b);
  }
  scale(f) {
    return new BABYLON.Color3(this.r * f, this.g * f, this.b * f);
  }
};
BABYLON.StandardMaterial ||= function () {
  this.diffuseColor = new BABYLON.Color3(1, 1, 1);
  this.emissiveColor = new BABYLON.Color3(0, 0, 0);
};
BABYLON.MeshBuilder ||= {
  CreateSphere: () => ({ thinInstanceEnablePicking: false, thinInstanceSetBuffer() {}, material: null }),
  CreateCylinder: () => ({ thinInstanceEnablePicking: false, thinInstanceSetBuffer() {}, material: null }),
};
BABYLON.PointerEventTypes ||= { POINTERDOWN: 1 };

function loadXYZ(name) {
  return fs.readFileSync(path.join(process.cwd(), 'public', 'molecules', name), 'utf-8');
}

function buildView(xyzText) {
  const parsed = parseXYZ(xyzText);
  const state = createMoleculeState();
  applyXYZToState(state, parsed);
  const bondSvc = createBondService(state);
  bondSvc.recomputeAndStore();
  const view = createMoleculeView({ onPointerObservable: { add() {} } }, state);
  return { state, view };
}

describe('x-visual-integration', () => {
  test('ROY atoms and bonds map to thin-instance groups', () => {
    const { state, view } = buildView(loadXYZ('roy.xyz'));
    const atomInstances = Array.from(view._internals.atomGroups.values()).reduce(
      (sum, group) => sum + group.indices.length,
      0
    );
    expect(atomInstances).toBe(state.positions.length);
    const bondInstances = Array.from(view._internals.bondGroups.values()).reduce(
      (sum, group) => sum + group.indices.length,
      0
    );
    expect(bondInstances).toBe(state.bonds.length);
  });

  test('Benzene atoms and bonds map to thin-instance groups', () => {
    const { state, view } = buildView(loadXYZ('benzene.xyz'));
    const atomInstances = Array.from(view._internals.atomGroups.values()).reduce(
      (sum, group) => sum + group.indices.length,
      0
    );
    expect(atomInstances).toBe(state.positions.length);
    const bondInstances = Array.from(view._internals.bondGroups.values()).reduce(
      (sum, group) => sum + group.indices.length,
      0
    );
    expect(bondInstances).toBe(state.bonds.length);
  });
});
