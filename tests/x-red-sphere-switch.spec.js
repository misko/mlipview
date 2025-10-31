import fs from 'fs';
import path from 'path';
import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.ts';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Minimal BABYLON scaffolding shared with other highlight-oriented x-tests.
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
  static Zero() {
    return new Vector3(0, 0, 0);
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
  static Dot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
  static Cross(a, b) {
    return new Vector3(
      a.y * b.z - a.z * b.y,
      a.z * b.x - a.x * b.z,
      a.x * b.y - a.y * b.x
    );
  }
}

class Quaternion {
  constructor(x = 0, y = 0, z = 0, w = 1) {
    this.x = x;
    this.y = y;
    this.z = z;
    this.w = w;
  }
  static Identity() {
    return new Quaternion();
  }
  static RotationAxis(axis, angle) {
    const half = angle / 2;
    const s = Math.sin(half);
    return new Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(half));
  }
}

class StandardMaterial {
  constructor() {
    this.diffuseColor = {};
    this.emissiveColor = {};
    this.specularColor = {};
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

function makeHighlightMesh(name) {
  return {
    name,
    position: Vector3.Zero(),
    scaling: new Vector3(1, 1, 1),
    rotationQuaternion: Quaternion.Identity(),
    thinInstanceSetBuffer() {},
    thinInstanceRefreshBoundingInfo() {},
    parent: null,
    isPickable: false,
    isVisible: false,
  };
}

const MeshBuilder = {
  CreateSphere: (name) => makeHighlightMesh(name),
  CreateCylinder: (name) => makeHighlightMesh(name),
  CreateLines: () => ({}),
};

Object.assign(global.BABYLON, {
  Vector3,
  Quaternion,
  StandardMaterial,
  Matrix,
  MeshBuilder,
  Axis: { X: new Vector3(1, 0, 0), Y: new Vector3(0, 1, 0), Z: new Vector3(0, 0, 1) },
});
global.BABYLON.Vector3.Up = () => new Vector3(0, 1, 0);

function makeScene() {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    getMeshByName() {
      return null;
    },
  };
}

function loadXYZ(name) {
  const file = path.join(process.cwd(), 'public', 'molecules', name);
  if (!fs.existsSync(file)) return null;
  return fs.readFileSync(file, 'utf8');
}

describe('red sphere artifact regression', () => {
  test('atom highlight hides after switching to smaller molecule', () => {
    const royText = loadXYZ('roy.xyz');
    const benzText = loadXYZ('benzene.xyz');
    if (!royText || !benzText) {
      console.warn('[x-red-sphere-switch] skipped (missing molecule fixtures)');
      return;
    }

    const royParsed = parseXYZ(royText);
    const state = createMoleculeState({
      elements: royParsed.elements,
      positions: royParsed.positions,
    });
    const scene = makeScene();
    const view = createMoleculeView(scene, state);
    const selection = createSelectionService(state);

    selection.clickAtom(0);
    const highlight = view?._internals?.highlight?.atom;
    expect(highlight).toBeTruthy();
    expect(highlight.isVisible).toBe(true);

    const benzParsed = parseXYZ(benzText);
    applyXYZToState(state, benzParsed);

    // Bond + position change events should force highlight to hide even though selection still references old index.
    expect(highlight.isVisible).toBe(false);
  });
});
