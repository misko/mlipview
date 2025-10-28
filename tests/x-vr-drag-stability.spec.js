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
  MeshBuilder: {
    CreateSphere: () => ({ thinInstanceEnablePicking: true, thinInstanceSetBuffer() {} }),
    CreateCylinder: () => ({ thinInstanceEnablePicking: true, thinInstanceSetBuffer() {} }),
  },
});

function makeScene(jitterFn) {
  return {
    meshes: [],
    onBeforeRenderObservable: { add() {} },
    getEngine: () => ({ getRenderingCanvas: () => null }),
    pickWithRay: () => ({ hit: false }),
    getControllerRay: jitterFn,
  };
}

describe('x-vr drag stability', () => {
  test('drag path remains smooth under jitter', () => {
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
    let step = 0;
    const jitter = () => {
      const dx = Math.sin(step * 3.7) * 0.02;
      const dy = Math.cos(step * 5.1) * 0.02;
      step++;
      return new Ray(new Vector3(0, 0, 0), new Vector3(dx, dy, -1).normalize(), 2);
    };
    const scene = makeScene(jitter);
    const vr = createVRSupport(scene, {
      picking: { molState: mol, selectionService: selection, manipulation },
    });
    vr.init();

    selection.clickAtom(1);
    const init = { ...mol.positions[1] };
    const planePoint = { ...init };
    const planeNormal = { x: 0, y: 1, z: 0 };
    const intersector = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
    manipulation.beginDrag(intersector, { planePoint, planeNormal });

    const path = [];
    for (let i = 0; i < 60; i++) {
      const ray = jitter();
      const point = {
        x: init.x + ray.direction.x * 0.3,
        y: init.y,
        z: init.z + ray.direction.z * 0.3,
      };
      const inter = () => point;
      manipulation.updateDrag(inter);
      path.push(mol.positions[1].z);
    }
    manipulation.endDrag();

    let maxDelta = 0;
    for (let i = 1; i < path.length; i++) {
      maxDelta = Math.max(maxDelta, Math.abs(path[i] - path[i - 1]));
    }
    expect(maxDelta).toBeLessThan(0.3);
  });
});
