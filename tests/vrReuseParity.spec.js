// Parity test: ensure initVRApp reuses existing viewer engine/scene and does not destroy atoms.
// We simulate a subset of the viewer environment: a mock engine/scene with meshes added.
// This is a structural test (no real Babylon rendering under jsdom).

function makeMockScene() {
  return {
    meshes: [],
    _renderCalls: 0,
    render() {
      this._renderCalls++;
    },
    createDefaultXRExperienceAsync: async () => ({ baseExperience: {} }),
  };
}
function makeMockEngine(scene) {
  return {
    _loops: 0,
    runRenderLoop(cb) {
      this._cb = cb;
    },
    step() {
      if (this._cb) {
        this._loops++;
        this._cb();
      }
    },
    resize() {},
    // minimal compatibility for code expecting engine reference
  };
}

describe('VR reuse parity', () => {
  test('initVRApp reuses existing global _viewer scene & engine', async () => {
    // Arrange a pseudo-viewer environment
    const scene = makeMockScene();
    // Pretend we have some atom meshes already
    for (let i = 0; i < 10; i++) scene.meshes.push({ name: 'atom_' + i });
    const engine = makeMockEngine(scene);
    global.window = global.window || {};
    window._viewer = { engine, scene };

    const { initVRApp } = await import('../public/vr/main-vr.js');
    const beforeMeshes = scene.meshes.length;
    const r = await initVRApp();
    expect(r.scene).toBe(scene);
    expect(r.engine).toBe(engine);
    expect(scene.meshes.length).toBe(beforeMeshes); // Should not add random meshes
  });
});
