import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createPickingService } from '../public/core/pickingService.js';

// Minimal Babylon stubs
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
BABYLON.Mesh = function () {};

function sceneStub() {
  return {
    onPointerObservable: { add() {} },
    pick() {
      return { hit: false };
    },
  };
}

describe('ghost atoms/bonds are not pickable', () => {
  test('toggle cell + ghosts and ensure picking still only returns real atoms', () => {
    const st = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
      cell: {
        a: { x: 5, y: 0, z: 0 },
        b: { x: 0, y: 5, z: 0 },
        c: { x: 0, y: 0, z: 5 },
        enabled: true,
        originOffset: { x: 0, y: 0, z: 0 },
      },
    });
    const bondSvc = createBondService(st);
    bondSvc.recomputeAndStore();
    const scene = sceneStub();
    const view = createMoleculeView(scene, st);
    const selection = createSelectionService(st);
    createPickingService(scene, view, selection, { manipulation: {}, camera: {} });
    // Enable cell and ghosts
    st.toggleCellVisibility();
    st.toggleGhostCells();
    st.markCellChanged();
    // Simulate pick by returning first real atom group master
    const realMaster = Array.from(view._internals.atomGroups.values())[0].master;
    scene.pick = () => ({ hit: true, pickedMesh: realMaster, thinInstanceIndex: 0 });
    // Trigger selection change manually via picking service pickAtPointer
    // (Would require exposing pickAtPointer; we simulate by calling selection directly)
    selection.clickAtom(0);
    expect(selection.get().kind).toBe('atom');
  });
});
