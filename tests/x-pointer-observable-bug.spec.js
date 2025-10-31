import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.ts';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createPickingService } from '../public/core/pickingService.js';

global.BABYLON = global.BABYLON || {};
const BABYLON = global.BABYLON;
BABYLON.PointerEventTypes = { POINTERDOWN: 1 };
BABYLON.Engine = function () {
  this.getRenderingCanvas = () => ({ getBoundingClientRect: () => ({ left: 0, top: 0 }) });
};
BABYLON.Scene = function () {
  this.onPointerObservable = {
    _listeners: [],
    add(fn) {
      this._listeners.push(fn);
    },
    notifyObservers(ev) {
      this._listeners.forEach((fn) => fn(ev));
    },
    notify(ev) {
      this.notifyObservers(ev);
    },
  };
};
BABYLON.Color4 = function () {};
BABYLON.Vector3 = function (x, y, z) {
  this.x = x;
  this.y = y;
  this.z = z;
};
BABYLON.Vector3.prototype.length = function () {
  return Math.hypot(this.x, this.y, this.z);
};
BABYLON.Vector3.Dot = (a, b) => a.x * b.x + a.y * b.y + a.z * b.z;
BABYLON.Vector3.Cross = (a, b) =>
  new BABYLON.Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
BABYLON.Vector3.prototype.normalize = function () {
  const L = this.length() || 1;
  this.x /= L;
  this.y /= L;
  this.z /= L;
  return this;
};
BABYLON.Vector3.prototype.normalizeToNew = function () {
  return new BABYLON.Vector3(this.x, this.y, this.z).normalize();
};
BABYLON.ArcRotateCamera = function () {
  return { attachControl() {} };
};
BABYLON.HemisphericLight = function () {};
BABYLON.StandardMaterial = function () {
  this.diffuseColor = { clone: () => ({}), scale: () => ({}) };
};
BABYLON.Color3 = function (r, g, b) {
  this.r = r;
  this.g = g;
  this.b = b;
  this.scale = (f) => new BABYLON.Color3(r * f, g * f, b * f);
  this.clone = () => new BABYLON.Color3(r, g, b);
};
BABYLON.MeshBuilder = {
  CreateSphere: () => ({ thinInstanceSetBuffer() {} }),
  CreateCylinder: () => ({ thinInstanceSetBuffer() {} }),
};
BABYLON.Matrix = class {
  constructor() {
    this.m = new Float32Array(16);
  }
  static Compose() {
    return new BABYLON.Matrix();
  }
};
BABYLON.Quaternion = { Identity: () => ({}), RotationAxis: () => ({}) };

test('x-pointer-observable handles notifyObservers-only implementation', async () => {
  const state = createMoleculeState({
    elements: ['H'],
    positions: [{ x: 0, y: 0, z: 0 }],
    bonds: [],
  });
  const canvas = { addEventListener() {}, getBoundingClientRect: () => ({ left: 0, top: 0 }) };
  const { createScene } = await import('../public/render/scene.js');
  const { scene } = await createScene(canvas);
  const view = createMoleculeView(scene, state);
  const selection = createSelectionService(state);
  scene.pick = () => ({
    hit: true,
    pickedMesh: Array.from(view._internals.atomGroups.values())[0].master,
    thinInstanceIndex: 0,
  });
  createPickingService(scene, view, selection);
  scene.onPointerObservable.notifyObservers({ type: BABYLON.PointerEventTypes.POINTERDOWN });
  expect(selection.get().kind).toBe('atom');
});

test('x-pointer-observable handles normalized notify alias', async () => {
  const state = createMoleculeState({
    elements: ['H', 'H'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
    ],
    bonds: [{ i: 0, j: 1 }],
  });
  const canvas = { addEventListener() {}, getBoundingClientRect: () => ({ left: 0, top: 0 }) };
  const { createScene } = await import('../public/render/scene.js');
  const { scene } = await createScene(canvas);
  const view = createMoleculeView(scene, state);
  const selection = createSelectionService(state);
  view.rebuildBonds();
  const firstBondGroup = Array.from(view._internals.bondGroups.values())[0];
  scene.pick = () => ({ hit: true, pickedMesh: firstBondGroup.master, thinInstanceIndex: 0 });
  createPickingService(scene, view, selection);
  scene.onPointerObservable.notify({ type: BABYLON.PointerEventTypes.POINTERDOWN });
  expect(selection.get().kind).toBe('bond');
});
