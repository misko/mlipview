import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createPickingService } from '../public/core/pickingService.js';

global.BABYLON = global.BABYLON || {};
const BABYLON = global.BABYLON;
BABYLON.PointerEventTypes ||= { POINTERDOWN: 1 };
BABYLON.StandardMaterial ||= function () {
  this.diffuseColor = { clone: () => ({}), scale: () => ({}) };
};
BABYLON.Color3 ||= function (r, g, b) {
  this.r = r;
  this.g = g;
  this.b = b;
  this.scale = (f) => new BABYLON.Color3(r * f, g * f, b * f);
  this.clone = () => new BABYLON.Color3(r, g, b);
};
BABYLON.MeshBuilder ||= {
  CreateSphere: () => ({ thinInstanceSetBuffer() {} }),
  CreateCylinder: () => ({ thinInstanceSetBuffer() {} }),
};
BABYLON.Matrix ||= class {
  constructor() {
    this.m = new Float32Array(16);
  }
  static Compose() {
    return new BABYLON.Matrix();
  }
};
BABYLON.Vector3 ||= class {
  constructor(x = 0, y = 0, z = 0) {
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
BABYLON.Quaternion ||= { Identity: () => ({}), RotationAxis: () => ({}) };

function makeScene(pickResultRef) {
  return {
    pointerX: 0,
    pointerY: 0,
    pick: () => pickResultRef.current,
    onPointerObservable: {
      add(fn) {
        this._fn = fn;
      },
      _fn: null,
      simulateDown() {
        this._fn?.({ type: BABYLON.PointerEventTypes.POINTERDOWN });
      },
    },
  };
}

describe('x-picking-integration', () => {
  test('pointer down over atom master selects atom', () => {
    const state = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const pickRef = { current: null };
    const scene = makeScene(pickRef);
    const view = createMoleculeView(scene, state);
    const selection = createSelectionService(state);
    createPickingService(scene, view, selection);
    const firstGroup = Array.from(view._internals.atomGroups.values())[0];
    scene.pick = () => ({ hit: true, pickedMesh: firstGroup.master, thinInstanceIndex: 0 });
    scene.onPointerObservable.simulateDown();
    expect(selection.get().kind).toBe('atom');
  });
});
