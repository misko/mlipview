export function installBabylonStub() {
  if (global.BABYLON && global.BABYLON.__isTestStub) {
    return global.BABYLON;
  }

  class Vector3 {
    constructor(x = 0, y = 0, z = 0) {
      this.x = x;
      this.y = y;
      this.z = z;
    }
    clone() {
      return new Vector3(this.x, this.y, this.z);
    }
    copyFrom(v) {
      this.x = v.x;
      this.y = v.y;
      this.z = v.z;
      return this;
    }
    add(v) {
      return new Vector3(this.x + v.x, this.y + v.y, this.z + v.z);
    }
    addInPlace(v) {
      this.x += v.x;
      this.y += v.y;
      this.z += v.z;
      return this;
    }
    subtract(v) {
      return new Vector3(this.x - v.x, this.y - v.y, this.z - v.z);
    }
    subtractInPlace(v) {
      this.x -= v.x;
      this.y -= v.y;
      this.z -= v.z;
      return this;
    }
    scale(f) {
      return new Vector3(this.x * f, this.y * f, this.z * f);
    }
    scaleInPlace(f) {
      this.x *= f;
      this.y *= f;
      this.z *= f;
      return this;
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
      return this.clone().normalize();
    }
    applyRotationQuaternionInPlace(q) {
      const x = this.x,
        y = this.y,
        z = this.z;
      const qx = q.x,
        qy = q.y,
        qz = q.z,
        qw = q.w;
      const ix = qw * x + qy * z - qz * y;
      const iy = qw * y + qz * x - qx * z;
      const iz = qw * z + qx * y - qy * x;
      const iw = -qx * x - qy * y - qz * z;
      this.x = ix * qw + iw * -qx + iy * -qz - iz * -qy;
      this.y = iy * qw + iw * -qy + iz * -qx - ix * -qz;
      this.z = iz * qw + iw * -qz + ix * -qy - iy * -qx;
      return this;
    }
    static Dot(a, b) {
      return a.x * b.x + a.y * b.y + a.z * b.z;
    }
    static Cross(a, b) {
      return new Vector3(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
    }
    static Zero() {
      return new Vector3(0, 0, 0);
    }
    static One() {
      return new Vector3(1, 1, 1);
    }
    static Up() {
      return new Vector3(0, 1, 0);
    }
    static Forward() {
      return new Vector3(0, 0, 1);
    }
    static Minimize(a, b) {
      return new Vector3(Math.min(a.x, b.x), Math.min(a.y, b.y), Math.min(a.z, b.z));
    }
    static Maximize(a, b) {
      return new Vector3(Math.max(a.x, b.x), Math.max(a.y, b.y), Math.max(a.z, b.z));
    }
  }

  class Quaternion {
    constructor(x = 0, y = 0, z = 0, w = 1) {
      this.x = x;
      this.y = y;
      this.z = z;
      this.w = w;
    }
    clone() {
      return new Quaternion(this.x, this.y, this.z, this.w);
    }
    copyFrom(q) {
      this.x = q.x;
      this.y = q.y;
      this.z = q.z;
      this.w = q.w;
      return this;
    }
    static Identity() {
      return new Quaternion();
    }
    static RotationAxis(axis, angle) {
      const s = Math.sin(angle / 2);
      return new Quaternion(axis.x * s, axis.y * s, axis.z * s, Math.cos(angle / 2));
    }
  }

  class Matrix {
    constructor() {
      this.m = new Float32Array(16);
    }
    static Identity() {
      const m = new Matrix();
      m.m[0] = m.m[5] = m.m[10] = m.m[15] = 1;
      return m;
    }
    static RotationAxis(axis, angle) {
      const m = Matrix.Identity();
      const s = Math.sin(angle),
        c = Math.cos(angle),
        t = 1 - c;
      const x = axis.x,
        y = axis.y,
        z = axis.z;
      m.m[0] = c + x * x * t;
      m.m[1] = x * y * t + z * s;
      m.m[2] = x * z * t - y * s;
      m.m[4] = y * x * t - z * s;
      m.m[5] = c + y * y * t;
      m.m[6] = y * z * t + x * s;
      m.m[8] = z * x * t + y * s;
      m.m[9] = z * y * t - x * s;
      m.m[10] = c + z * z * t;
      return m;
    }
    static Compose(scale, rot, trans) {
      const mat = new Matrix();
      mat.m[0] = scale.x;
      mat.m[5] = scale.y;
      mat.m[10] = scale.z;
      mat.m[12] = trans.x;
      mat.m[13] = trans.y;
      mat.m[14] = trans.z;
      mat.m[15] = 1;
      if (rot) {
        // Apply very lightweight quaternion-to-matrix (not exact but adequate for tests)
        const { x, y, z, w } = rot;
        const xx = x * x,
          yy = y * y,
          zz = z * z;
        const xy = x * y,
          xz = x * z,
          yz = y * z;
        const wx = w * x,
          wy = w * y,
          wz = w * z;
        mat.m[0] = 1 - 2 * (yy + zz);
        mat.m[1] = 2 * (xy + wz);
        mat.m[2] = 2 * (xz - wy);
        mat.m[4] = 2 * (xy - wz);
        mat.m[5] = 1 - 2 * (xx + zz);
        mat.m[6] = 2 * (yz + wx);
        mat.m[8] = 2 * (xz + wy);
        mat.m[9] = 2 * (yz - wx);
        mat.m[10] = 1 - 2 * (xx + yy);
      }
      return mat;
    }
  }

  class Color3 {
    constructor(r = 0, g = 0, b = 0) {
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

  class Color4 {
    constructor(r = 0, g = 0, b = 0, a = 1) {
      this.r = r;
      this.g = g;
      this.b = b;
      this.a = a;
    }
    clone() {
      return new Color4(this.r, this.g, this.b, this.a);
    }
  }

  class StandardMaterial {
    constructor(name) {
      this.name = name;
      this.disableLighting = false;
      this.emissiveColor = new Color3(0, 0, 0);
      this.alpha = 1;
    }
  }

  class Ray {
    constructor(origin = Vector3.Zero(), direction = Vector3.Forward(), length = 1) {
      this.origin = origin;
      this.direction = direction;
      this.length = length;
    }
  }

  const MeshBuilder = {
    CreateSphere(name, _opts, scene) {
      const mesh = {
        name,
        position: Vector3.Zero(),
        scaling: Vector3.One(),
        thinInstanceEnablePicking: true,
        thinInstanceSetBuffer() {},
        isPickable: true,
        dispose() {},
      };
      if (scene?.meshes) scene.meshes.push(mesh);
      return mesh;
    },
    CreateCylinder(name, _opts, scene) {
      const mesh = {
        name,
        position: Vector3.Zero(),
        scaling: Vector3.One(),
        thinInstanceEnablePicking: true,
        thinInstanceSetBuffer() {},
        isPickable: true,
        parent: null,
        material: null,
        dispose() {},
      };
      if (scene?.meshes) scene.meshes.push(mesh);
      return mesh;
    },
    CreateLines(name, _opts, scene) {
      const mesh = {
        name,
        color: new Color3(1, 1, 1),
        alpha: 1,
        dispose() {},
      };
      if (scene?.meshes) scene.meshes.push(mesh);
      return mesh;
    },
  };

  class TransformNode {
    constructor(name) {
      this.name = name;
      this.position = Vector3.Zero();
      this.scaling = Vector3.One();
      this.rotationQuaternion = Quaternion.Identity();
      this._children = [];
    }
    getClassName() {
      return 'TransformNode';
    }
    addChild(node) {
      node.parent = this;
      this._children.push(node);
    }
    getChildren() {
      return this._children.slice();
    }
    getChildMeshes() {
      return this._children.filter((c) => c.getClassName && c.getClassName() === 'Mesh');
    }
  }

  class Engine {
    constructor(canvas) {
      this.canvas = canvas;
      this._loop = null;
    }
    setHardwareScalingLevel() {}
    getRenderingCanvas() {
      return this.canvas;
    }
    runRenderLoop(cb) {
      this._loop = cb;
    }
    step() {
      if (this._loop) this._loop();
    }
    resize() {}
  }

  class Scene {
    constructor(engine) {
      this.engine = engine;
      this.meshes = [];
      this.transformNodes = [];
      this.clearColor = new Color4();
    }
    getEngine() {
      return this.engine;
    }
    render() {}
    createDefaultXRExperienceAsync() {
      return Promise.resolve({
        input: { controllers: [], onControllerAddedObservable: { add() {} } },
        baseExperience: {},
      });
    }
  }

  class ArcRotateCamera {
    constructor(name, alpha, beta, radius, target) {
      this.name = name;
      this.alpha = alpha;
      this.beta = beta;
      this.radius = radius;
      this.target = target;
      this.position = Vector3.Zero();
      this.fov = Math.PI / 4;
    }
    attachControl() {}
    getDirection() {
      return Vector3.Forward();
    }
  }

  const stub = {
    Vector3,
    Quaternion,
    Matrix,
    Color3,
    Color4,
    StandardMaterial,
    MeshBuilder,
    TransformNode,
    ArcRotateCamera,
    Engine,
    Scene,
    Ray,
    Axis: {
      X: new Vector3(1, 0, 0),
      Y: new Vector3(0, 1, 0),
      Z: new Vector3(0, 0, 1),
    },
    WebXRFeatureName: {
      POINTER_SELECTION: 'pointer-selection',
    },
    GUI: {},
    __isTestStub: true,
  };

  global.BABYLON = stub;
  return stub;
}
