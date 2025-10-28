// Minimal Babylon-like math/mesh stub for Jest environments.
// Provides just enough implementation for moleculeView force and highlight tests.

function createVector3Class() {
  return class Vector3 {
    constructor(x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    length() {
      return Math.hypot(this.x, this.y, this.z);
    }
    addInPlace(v) {
      this.x += v.x;
      this.y += v.y;
      this.z += v.z;
      return this;
    }
    subtractInPlace(v) {
      this.x -= v.x;
      this.y -= v.y;
      this.z -= v.z;
      return this;
    }
    scaleInPlace(f) {
      this.x *= f;
      this.y *= f;
      this.z *= f;
      return this;
    }
    clone() {
      return new Vector3(this.x, this.y, this.z);
    }
    normalize() {
      const L = this.length() || 1;
      this.x /= L;
      this.y /= L;
      this.z /= L;
      return this;
    }
    normalizeToNew() {
      return this.clone().normalize();
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
    static TransformCoordinates(v, matrix) {
      if (!matrix || !matrix.m) return v.clone();
      const m = matrix.m;
      const x = v.x,
        y = v.y,
        z = v.z;
      return new Vector3(
        m[0] * x + m[4] * y + m[8] * z + m[12],
        m[1] * x + m[5] * y + m[9] * z + m[13],
        m[2] * x + m[6] * y + m[10] * z + m[14]
      );
    }
    static Zero() {
      return new Vector3(0, 0, 0);
    }
  };
}

function createQuaternionClass() {
  return class Quaternion {
    constructor(x = 0, y = 0, z = 0, w = 1) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.w = w;
    }
    static Identity() {
      return new Quaternion(0, 0, 0, 1);
    }
    static RotationAxis(axis, angle) {
      const half = angle / 2;
      const s = Math.sin(half);
      const norm = Math.hypot(axis?.x ?? 0, axis?.y ?? 0, axis?.z ?? 0) || 1;
      return new Quaternion(
        (axis?.x ?? 0) / norm * s,
        (axis?.y ?? 0) / norm * s,
        (axis?.z ?? 0) / norm * s,
        Math.cos(half)
      );
    }
  };
}

function createMatrixClass(Vector3Class, QuaternionClass) {
  return class Matrix {
    constructor() {
      this.m = new Float32Array(16);
      Matrix.IdentityToRef(this);
    }
    static Identity() {
      const mat = new Matrix();
      Matrix.IdentityToRef(mat);
      return mat;
    }
    static IdentityToRef(target) {
      const m = target.m || (target.m = new Float32Array(16));
      m.fill(0);
      m[0] = m[5] = m[10] = m[15] = 1;
    }
    static Compose(scale, rotation, translation) {
      const mat = new Matrix();
      // Very small subset: scale diagonal + translation; ignore rotation for tests.
      const sx = scale?.x ?? scale?.y ?? scale?.z ?? 1;
      const sy = scale?.y ?? sx;
      const sz = scale?.z ?? sx;
      mat.m[0] = sx;
      mat.m[5] = sy;
      mat.m[10] = sz;
      mat.m[12] = translation?.x ?? 0;
      mat.m[13] = translation?.y ?? 0;
      mat.m[14] = translation?.z ?? 0;
      mat.__scale = new Vector3Class(sx, sy, sz);
      mat.__rotation = rotation instanceof QuaternionClass ? rotation : QuaternionClass.Identity();
      mat.__translation = translation instanceof Vector3Class ? translation.clone() : new Vector3Class(
        translation?.x ?? 0,
        translation?.y ?? 0,
        translation?.z ?? 0
      );
      return mat;
    }
  };
}

class Color3 {
  constructor(r = 1, g = 0, b = 0) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  clone() {
    return new Color3(this.r, this.g, this.b);
  }
  scale(f) {
    return new Color3(this.r * f, this.g * f, this.b * f);
  }
}

class StandardMaterial {
  constructor() {
    this.diffuseColor = new Color3();
    this.emissiveColor = new Color3();
    this.specularColor = new Color3();
    this.alpha = 1;
    this.disableLighting = false;
  }
}

class Mesh {
  constructor() {
    this._buffers = {};
    this.material = null;
    this.isVisible = true;
    this.thinInstanceEnablePicking = false;
  }
  thinInstanceSetBuffer(kind, arr) {
    this._buffers[kind] = arr;
  }
  setEnabled(on) {
    this.isVisible = !!on;
  }
  dispose() {}
}

const MeshBuilder = {
  CreateSphere: () => new Mesh(),
  CreateCylinder: () => new Mesh(),
};

export function ensureBabylonStub() {
  if (global.BABYLON && global.BABYLON.__isTestStub) return global.BABYLON;

  const Vector3 = createVector3Class();
  const Quaternion = createQuaternionClass();
  const Matrix = createMatrixClass(Vector3, Quaternion);

  global.BABYLON = {
    __isTestStub: true,
    Vector3,
    Quaternion,
    Matrix,
    Color3,
    StandardMaterial,
    MeshBuilder,
    TransformNode: class {},
    Material: { MATERIAL_ALPHABLEND: 2 },
  };

  return global.BABYLON;
}

export default { ensureBabylonStub };
