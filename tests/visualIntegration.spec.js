// Visual integration style test using stubs: ensures atom and bond thin instance bookkeeping
// matches the molecule state after loading ROY and Benzene.
import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import fs from 'fs';
import path from 'path';

// Minimal Babylon stub for matrices & materials sufficient for moleculeView logic.
global.BABYLON = (function () {
  class Vector3 {
    constructor(x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    static Dot(a, b) {
      return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    static Cross(a, b) {
      return new Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
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
  }
  class Quaternion {
    static Identity() {
      return new Quaternion();
    }
    static RotationAxis() {
      return new Quaternion();
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
  class Color3 {
    constructor(r, g, b) {
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
      this.emissiveColor = new Color3(0, 0, 0);
    }
  }
  class MeshBuilder {
    static CreateSphere(name) {
      return { name, thinInstanceEnablePicking: false, thinInstanceSetBuffer() {}, material: null };
    }
    static CreateCylinder(name) {
      return { name, thinInstanceEnablePicking: false, thinInstanceSetBuffer() {}, material: null };
    }
  }
  const PointerEventTypes = { POINTERDOWN: 1 };
  return { Vector3, Quaternion, Matrix, Color3, StandardMaterial, MeshBuilder, PointerEventTypes };
})();

function loadXYZFile(filename) {
  return fs.readFileSync(path.join(process.cwd(), 'public', 'molecules', filename), 'utf-8');
}

function buildViewFromXYZ(xyzText) {
  const parsed = parseXYZ(xyzText);
  const state = createMoleculeState();
  applyXYZToState(state, parsed);
  // compute bonds
  const bondSvc = createBondService(state);
  bondSvc.recomputeAndStore();
  const sceneStub = { onPointerObservable: { add() {} } }; // minimal scene stub
  const view = createMoleculeView(sceneStub, state);
  return { state, view };
}

describe('visual integration atoms & bonds', () => {
  test('ROY atoms & bonds represented', () => {
    const { state, view } = buildViewFromXYZ(loadXYZFile('roy.xyz'));
    const totalAtomsInGroups = Array.from(view._internals.atomGroups.values()).reduce(
      (s, g) => s + g.indices.length,
      0
    );
    expect(totalAtomsInGroups).toBe(state.positions.length);
    const totalBondsInGroups = Array.from(view._internals.bondGroups.values()).reduce(
      (s, g) => s + g.indices.length,
      0
    );
    expect(totalBondsInGroups).toBe(state.bonds.length);
  });
  test('Benzene atoms & bonds represented', () => {
    const { state, view } = buildViewFromXYZ(loadXYZFile('benzene.xyz'));
    const totalAtomsInGroups = Array.from(view._internals.atomGroups.values()).reduce(
      (s, g) => s + g.indices.length,
      0
    );
    expect(totalAtomsInGroups).toBe(state.positions.length);
    const totalBondsInGroups = Array.from(view._internals.bondGroups.values()).reduce(
      (s, g) => s + g.indices.length,
      0
    );
    expect(totalBondsInGroups).toBe(state.bonds.length);
  });
});
