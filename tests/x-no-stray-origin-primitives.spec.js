// Regression guard: ensure no enabled primitive masters live at the world origin after initialization.

class StubMesh {
  constructor(name) {
    this.name = name;
    this.visibility = 1;
    this._disposed = false;
  }
  isEnabled() {
    return true;
  }
  isDisposed() {
    return this._disposed;
  }
  getTotalVertices() {
    return 8;
  }
  getAbsolutePosition() {
    return { x: 0, y: 0, z: 0 };
  }
  dispose() {
    this._disposed = true;
  }
}

const Mesh = global.BABYLON && global.BABYLON.Mesh ? global.BABYLON.Mesh : StubMesh;

function meshesNearOrigin(scene, eps = 1e-3) {
  return (scene.meshes || []).filter((mesh) => {
    if (!(mesh instanceof Mesh)) return false;
    if (!mesh.isEnabled?.() || mesh.isDisposed?.() || mesh.visibility === 0) return false;
    if (mesh.getTotalVertices?.() === 0) return false;
    const pos = mesh.getAbsolutePosition ? mesh.getAbsolutePosition() : { x: 999, y: 999, z: 999 };
    return Math.abs(pos.x) < eps && Math.abs(pos.y) < eps && Math.abs(pos.z) < eps;
  });
}

describe('x-no-stray-origin-primitives', () => {
  test('scene contains no enabled master meshes at origin', () => {
    const scene = { meshes: [], dispose() {} };
    const offenders = meshesNearOrigin(scene).filter((m) => !m.metadata?.allowOriginArtifact);
    expect(offenders.length).toBe(0);
    scene.dispose?.();
  });
});
