import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createManipulationService } from '../public/domain/manipulationService.js';
import { createVRSupport } from '../public/vr/setup.js';

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
  subtract(v) {
    return new Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
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
  add(v) {
    return new Vector3(this.x + v.x, this.y + v.y, this.z + v.z);
  }
  scale(f) {
    return new Vector3(this.x * f, this.y * f, this.z * f);
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
}

class Ray {
  constructor(origin, direction, length = 100) {
    this.origin = origin;
    this.direction = direction;
    this.length = length;
  }
}

Object.assign(global.BABYLON, {
  Vector3,
  Quaternion,
  Ray,
  Axis: { X: new Vector3(1, 0, 0), Y: new Vector3(0, 1, 0), Z: new Vector3(0, 0, 1) },
  Matrix: class {
    constructor() {
      this.m = new Float32Array(16);
    }
    static Identity() {
      return new global.BABYLON.Matrix();
    }
  },
  MeshBuilder: {
    CreateSphere: () => ({ isPickable: true, thinInstanceEnablePicking: true, thinInstanceSetBuffer() {} }),
    CreateCylinder: () => ({ isPickable: true, thinInstanceEnablePicking: true, thinInstanceSetBuffer() {} }),
  },
});

function makeScene() {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    pickWithRay: () => ({ hit: false }),
    getEngine: () => ({ getRenderingCanvas: () => null }),
  };
}

describe('x-vr atom drag initiation', () => {
  beforeAll(() => {
    jest.useFakeTimers();
  });
  afterAll(() => {
    jest.useRealTimers();
  });

  test('selection + hold begins drag', () => {
    const mol = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [],
    });
    const selection = createSelectionService(mol);
    const manipulation = createManipulationService(mol, {});
    const scene = makeScene();
    const vr = createVRSupport(scene, {
      picking: { molState: mol, selectionService: selection, manipulation },
    });
    vr.init();

    selection.clickAtom(1);
    expect(mol.selection.kind).toBe('atom');

    const intersector = () => mol.positions[1];
    const started = manipulation.beginDrag(intersector, { planePoint: mol.positions[1], planeNormal: { x: 0, y: 1, z: 0 } });
    expect(started).toBe(true);

    const pre = { ...mol.positions[1] };
    const offset = { x: 0.2, y: 0, z: 0.3 };
    const moveIntersector = () => ({ x: pre.x + offset.x, y: pre.y + offset.y, z: pre.z + offset.z });
    manipulation.updateDrag(moveIntersector);
    const moved = mol.positions[1];
    expect(moved.x).toBeCloseTo(pre.x + offset.x, 6);
    expect(moved.z).toBeCloseTo(pre.z + offset.z, 6);

    manipulation.endDrag();
  });
});
