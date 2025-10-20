import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createPickingService } from '../public/core/pickingService.js';

// Stub minimal Babylon objects required for picking logic
global.BABYLON = global.BABYLON || {};
if (!BABYLON.PointerEventTypes) BABYLON.PointerEventTypes = { POINTERDOWN: 1 };
if (!BABYLON.StandardMaterial)
  BABYLON.StandardMaterial = function () {
    this.diffuseColor = { clone: () => ({}), scale: () => ({}) };
  };
if (!BABYLON.Color3)
  BABYLON.Color3 = function (r, g, b) {
    this.r = r;
    this.g = g;
    this.b = b;
    this.scale = (f) => new BABYLON.Color3(r * f, g * f, b * f);
    this.clone = () => new BABYLON.Color3(r, g, b);
  };
if (!BABYLON.MeshBuilder)
  BABYLON.MeshBuilder = {
    CreateSphere: () => ({ thinInstanceSetBuffer() {} }),
    CreateCylinder: () => ({ thinInstanceSetBuffer() {} }),
  };
if (!BABYLON.Matrix)
  BABYLON.Matrix = class {
    constructor() {
      this.m = new Float32Array(16);
    }
    static Compose() {
      return new BABYLON.Matrix();
    }
  };
if (!BABYLON.Vector3)
  BABYLON.Vector3 = class {
    constructor(x, y, z) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    length() {
      return Math.hypot(this.x, this.y, this.z);
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
    normalizeToNew() {
      const L = this.length() || 1;
      return new BABYLON.Vector3(this.x / L, this.y / L, this.z / L);
    }
  };
if (!BABYLON.Quaternion) BABYLON.Quaternion = { Identity: () => ({}), RotationAxis: () => ({}) };

function makeSceneStub(pickResultRef) {
  const scene = {
    pointerX: 0,
    pointerY: 0,
    pick: () => pickResultRef.current,
    onPointerObservable: {
      add(fn) {
        this._fn = fn;
      },
      _fn: null,
      simulateDown() {
        this._fn && this._fn({ type: BABYLON.PointerEventTypes.POINTERDOWN });
      },
    },
  };
  return scene;
}

describe('picking integration (selection + view)', () => {
  test('selects atom on pointer down', () => {
    const state = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const scene = makeSceneStub({ current: null });
    const view = createMoleculeView(scene, state);
    const selection = createSelectionService(state);
    const picking = createPickingService(scene, view, selection);
    // fabricate a pick pointing to first atom group master
    const firstGroup = Array.from(view._internals.atomGroups.values())[0];
    scene.pick = () => ({ hit: true, pickedMesh: firstGroup.master, thinInstanceIndex: 0 });
    scene.onPointerObservable.simulateDown();
    expect(selection.get().kind).toBe('atom');
  });
});
