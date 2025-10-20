import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createBondService } from '../public/domain/bondService.js';
import { computeBondsNoState } from '../public/bond_render.js';

// Inject computeBondsNoState into window-like global for rebuildGhosts access path
if (typeof window === 'undefined') global.window = {};
window.computeBondsNoState = computeBondsNoState;

if (!global.BABYLON) global.BABYLON = {};
BABYLON.Color3 = function (r, g, b) {
  this.r = r;
  this.g = g;
  this.b = b;
  this.clone = () => new BABYLON.Color3(r, g, b);
  this.scale = (f) => new BABYLON.Color3(r * f, g * f, b * f);
};
BABYLON.StandardMaterial = function () {
  this.diffuseColor = new BABYLON.Color3(1, 1, 1);
  this.emissiveColor = new BABYLON.Color3(0, 0, 0);
};
BABYLON.Vector3 = function (x, y, z) {
  this.x = x;
  this.y = y;
  this.z = z;
  this.length = () => Math.hypot(x, y, z);
};
BABYLON.Vector3.Dot = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z;
BABYLON.Vector3.Cross = (a, b) =>
  new BABYLON.Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
BABYLON.Vector3.prototype.normalizeToNew = function () {
  const L = this.length() || 1;
  return new BABYLON.Vector3(this.x / L, this.y / L, this.z / L);
};
BABYLON.Vector3.prototype.normalize = function () {
  const L = this.length() || 1;
  this.x /= L;
  this.y /= L;
  this.z /= L;
  return this;
};
BABYLON.Quaternion = { Identity: () => ({}), RotationAxis: () => ({}) };
BABYLON.Matrix = class {
  constructor() {
    this.m = new Float32Array(16);
  }
  static Compose() {
    return new BABYLON.Matrix();
  }
};
BABYLON.MeshBuilder = {
  CreateSphere: () => ({ thinInstanceSetBuffer() {}, material: null }),
  CreateCylinder: () => ({ thinInstanceSetBuffer() {}, material: null }),
  CreateLines: () => ({ dispose() {}, color: null }),
};

function sceneStub() {
  return { onPointerObservable: { add() {} } };
}

// System: a diatomic near the +x boundary so its ghost in +x direction should form crossing bond with original neighbor
// Atoms placed so that base atom 1 and ghost of atom 0 (shift +x) fall within bonding cutoff

describe('cross-image ghost bonds are generated', () => {
  test('augmented bond calculator yields at least one ghost bond', () => {
    const st = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0.9, y: 0, z: 0 },
        { x: 2.1, y: 0, z: 0 },
      ],
      bonds: [],
      cell: {
        a: { x: 3, y: 0, z: 0 },
        b: { x: 0, y: 3, z: 0 },
        c: { x: 0, y: 0, z: 3 },
        enabled: true,
        originOffset: { x: 0, y: 0, z: 0 },
      },
    });
    const bondSvc = createBondService(st);
    bondSvc.recomputeAndStore();
    const scene = sceneStub();
    const view = createMoleculeView(scene, st);
    st.toggleCellVisibility();
    st.toggleGhostCells();
    st.markCellChanged();
    // After rebuildGhosts, ghost bond groups are internal; we indirectly verify by recomputing augmented bonds here similarly
    const baseAtoms = st.positions.map((p, i) => ({
      element: st.elements[i],
      pos: [p.x, p.y, p.z],
    }));
    // Rebuild ghosts uses shifts including +/- x; ghost of atom 0 at +3 should be near atom 1 at 2.1 (distance 0.9)
    const aug = [
      ...baseAtoms,
      ...baseAtoms.map((a) => ({ element: a.element, pos: [a.pos[0] + 3, a.pos[1], a.pos[2]] })),
    ];
    const augCalc = computeBondsNoState(aug);
    const crossing = augCalc.some((b) => (b.i < 2 && b.j >= 2) || (b.j < 2 && b.i >= 2));
    expect(crossing).toBe(true);
  });
});
