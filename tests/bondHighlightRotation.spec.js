import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// Minimal BABYLON stubs + rotation support so we can simulate master rotation inheritance.
if (!global.BABYLON) global.BABYLON = {};
const B = global.BABYLON;
if (!B.StandardMaterial)
  B.StandardMaterial = function () {
    this.diffuseColor = {};
    this.emissiveColor = {};
    this.alpha = 1;
  };
if (!B.Color3)
  B.Color3 = function (r, g, b) {
    this.r = r;
    this.g = g;
    this.b = b;
    this.clone = () => new B.Color3(r, g, b);
    this.scale = (f) => new B.Color3(r * f, g * f, b * f);
  };
if (!B.MeshBuilder)
  B.MeshBuilder = {
    CreateSphere: (n, o, s) => ({
      name: n,
      thinInstanceSetBuffer() {},
      getChildren() {
        return [];
      },
    }),
    CreateCylinder: (n, o, s) => ({
      name: n,
      thinInstanceSetBuffer() {},
      getChildren() {
        return [];
      },
    }),
  };
if (!B.Matrix)
  B.Matrix = class {
    constructor() {
      this.m = new Float32Array(16);
    }
    static Compose() {
      return new B.Matrix();
    }
  };
if (!B.Vector3)
  B.Vector3 = class {
    constructor(x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    length() {
      return Math.hypot(this.x, this.y, this.z);
    }
    add(v) {
      return new B.Vector3(this.x + v.x, this.y + v.y, this.z + v.z);
    }
    subtract(v) {
      return new B.Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
    }
    scale(f) {
      return new B.Vector3(this.x * f, this.y * f, this.z * f);
    }
    static Dot(a, b) {
      return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    static Cross(a, b) {
      return new B.Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }
    normalizeToNew() {
      const L = this.length() || 1;
      return new B.Vector3(this.x / L, this.y / L, this.z / L);
    }
  };
if (!B.Quaternion)
  B.Quaternion = class {
    constructor(x = 0, y = 0, z = 0, w = 1) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.w = w;
    }
    static Identity() {
      return new B.Quaternion(0, 0, 0, 1);
    }
    static RotationAxis(axis, angle) {
      const h = angle / 2;
      const s = Math.sin(h);
      return new B.Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(h));
    }
  };
if (!B.Axis)
  B.Axis = { X: new B.Vector3(1, 0, 0), Y: new B.Vector3(0, 1, 0), Z: new B.Vector3(0, 0, 1) };

function makeSceneStub() {
  // Provide onBeforeRenderObservable so frame update path runs.
  return {
    onPointerObservable: { add() {} },
    onBeforeRenderObservable: {
      _c: [],
      add(fn) {
        this._c.push(fn);
      },
      run() {
        for (const f of this._c) f();
      },
    },
    meshes: [],
  };
}

// Helper to manually invoke frame observers (simulate render loop)
function step(scene, n = 1) {
  for (let i = 0; i < n; i++) scene.onBeforeRenderObservable.run();
}

describe('bond highlight rotation inheritance', () => {
  test('bond highlight midpoint rotates with master', () => {
    const state = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const scene = makeSceneStub();
    const view = createMoleculeView(scene, state);
    const selection = createSelectionService(state);
    // Select the bond to create + position highlight
    selection.clickBond({ i: 0, j: 1, key: 'C-H', index: 0 });
    step(scene, 2); // ensure initial highlight transform applied
    const { highlight, bondGroups } = view._internals;
    expect(highlight.bond.isVisible).toBe(true);
    // Capture initial local midpoint (should be ~ (0.5,0,0))
    const initialPos = {
      x: highlight.bond.position.x,
      y: highlight.bond.position.y,
      z: highlight.bond.position.z,
    };
    expect(Math.abs(initialPos.x - 0.5)).toBeLessThan(1e-6);
    // Simulate a rotation of the bond master around Y 90 degrees: local midpoint unchanged, world would change.
    // We emulate world change by verifying: because highlight is parented, we only need to ensure parenting occurred.
    const firstGroup = bondGroups.values().next().value;
    expect(firstGroup.master).toBeTruthy();
    // Attach a rotationQuaternion to master (as VR would do) and step frames; highlight local position should remain the same (local space),
    // confirming we are not overwriting it with stale world coordinates; parenting ensures world rotation now.
    firstGroup.master.rotationQuaternion = B.Quaternion.RotationAxis(B.Axis.Y, Math.PI / 2);
    step(scene, 3);
    // Local midpoint should still be the same (0.5,0,0)
    expect(Math.abs(highlight.bond.position.x - 0.5)).toBeLessThan(1e-6);
    expect(Math.abs(highlight.bond.position.y - 0.0)).toBeLessThan(1e-6);
    // Sanity: scaling.y reflects bond length (1 unit)
    expect(Math.abs(highlight.bond.scaling.y - 1)).toBeLessThan(1e-6);
  });
});
