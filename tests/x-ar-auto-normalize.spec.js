import { computeMoleculeWorldBounds } from '../public/vr/vr-utils.js';

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
  static Minimize(a, b) {
    return new Vector3(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.min(a.z, b.z));
  }
  static Maximize(a, b) {
    return new Vector3(Math.max(a.x, b.x), Math.max(a.y, b.y), Math.max(a.z, b.z));
  }
  add(v) {
    return new Vector3(this.x + v.x, this.y + v.y, this.z + v.z);
  }
  subtract(v) {
    return new Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
  }
  scale(k) {
    return new Vector3(this.x * k, this.y * k, this.z * k);
  }
  length() {
    return Math.hypot(this.x, this.y, this.z);
  }
}

class Mesh {
  constructor(name, min, max) {
    this.name = name;
    this._boundingInfo = {
      boundingBox: {
        minimumWorld: min,
        maximumWorld: max,
      },
    };
    this.parent = null;
  }
  refreshBoundingInfo() {}
  getBoundingInfo() {
    return this._boundingInfo;
  }
}

class TransformNode {
  constructor(name) {
    this.name = name;
    this.position = new Vector3();
    this.scaling = new Vector3(1, 1, 1);
    this._children = [];
  }
  addChild(mesh) {
    mesh.parent = this;
    this._children.push(mesh);
  }
  getChildMeshes() {
    return this._children;
  }
}

Object.assign(global.BABYLON, { Vector3, TransformNode });

describe('x-ar-auto-normalize', () => {
  test('bounds computed from root children when only molecule_root is present', () => {
    const scene = { meshes: [], transformNodes: [] };
    const root = new global.BABYLON.TransformNode('molecule_root');
    scene.transformNodes.push(root);
    const m1 = new Mesh('atom_O', new Vector3(-1, -1, -1), new Vector3(0, 0, 0));
    const m2 = new Mesh('atom_H', new Vector3(0, 0, 0), new Vector3(1, 1, 1));
    root.addChild(m1);
    root.addChild(m2);
    scene.meshes.push(m1, m2);

    const bounds = computeMoleculeWorldBounds(scene, root);
    expect(bounds).toBeTruthy();
    expect(bounds.min.x).toBe(-1);
    expect(bounds.max.z).toBe(1);
  });

  test('returns null when geometry is absent', () => {
    const scene = { meshes: [], transformNodes: [] };
    const root = new global.BABYLON.TransformNode('molecule_root');
    scene.transformNodes.push(root);

    const result = computeMoleculeWorldBounds(scene, root);
    expect(result).toBeNull();
  });
});

