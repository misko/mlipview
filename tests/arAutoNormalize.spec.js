// AR auto-normalization regression test
// Ensures molecule is placed forward in AR even when only molecule_root is present as master.

import { computeMoleculeWorldBounds } from '../public/vr/vr-utils.js';

// Babylon stubs (minimal for bounds math)
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
  scale(f) {
    return new Vector3(this.x * f, this.y * f, this.z * f);
  }
  subtract(v) {
    return new Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
  }
  length() {
    return Math.hypot(this.x, this.y, this.z);
  }
}
class Quaternion {
  constructor(x = 0, y = 0, z = 0, w = 1) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
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
  getChildMeshes() {
    return this._children;
  }
  addChild(m) {
    m.parent = this;
    this._children.push(m);
  }
}
Object.assign(BABYLON, { Vector3, Quaternion, TransformNode });

describe('ar auto-normalization', () => {
  test('computes bounds from root children if masters is only root', () => {
    const scene = { meshes: [], transformNodes: [] };
    const root = new BABYLON.TransformNode('molecule_root');
    scene.transformNodes.push(root);
    // Add two child meshes with bounding info
    const m1 = new Mesh('atom_O', new Vector3(-1, -1, -1), new Vector3(0, 0, 0));
    const m2 = new Mesh('atom_H', new Vector3(0, 0, 0), new Vector3(1, 1, 1));
    root.addChild(m1);
    root.addChild(m2);
    scene.meshes.push(m1, m2);
    // Should compute min (-1,-1,-1), max (1,1,1)
    const { min, max } = computeMoleculeWorldBounds(scene, root);
    expect(min.x).toBe(-1);
    expect(max.z).toBe(1);
  });
  test('returns null if no geometry', () => {
    const scene = { meshes: [], transformNodes: [] };
    const root = new BABYLON.TransformNode('molecule_root');
    scene.transformNodes.push(root);
    // No children
    const result = computeMoleculeWorldBounds(scene, root);
    expect(result).toBe(null);
  });
});
