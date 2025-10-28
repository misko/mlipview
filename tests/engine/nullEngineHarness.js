// NullEngine harness for Babylon scene graph logic tests (no real rendering)
// Intentionally lightweight: if full @babylonjs/core not present, this file can be adapted.

let NullEngine, Scene, MeshBuilder;
try {
  // Attempt modular imports (user may later add @babylonjs/core)
  ({ NullEngine } = require('@babylonjs/core/Engines/nullEngine'));
  ({ Scene } = require('@babylonjs/core/scene'));
  ({ MeshBuilder } = require('@babylonjs/core/Meshes/meshBuilder')); // side-effect import pattern
} catch (e) {
  // Fallback: use existing global BABYLON mocks if core not installed.
  if (global.BABYLON) {
    NullEngine = class {
      constructor() {
        this.fps = 60;
      }
      runRenderLoop(cb) {
        cb();
      }
      stopRenderLoop() {}
      getFps() {
        return this.fps;
      }
      dispose() {}
    };
    Scene = class {
      constructor() {
        this.meshes = [];
      }
      render() {}
      dispose() {}
    };
    MeshBuilder = global.BABYLON.MeshBuilder;
  } else {
    throw new Error(
      'NullEngine harness requires either @babylonjs/core or the global BABYLON mock.'
    );
  }
}

function createNullScene() {
  const engine = new NullEngine();
  const scene = new Scene(engine);
  return { engine, scene };
}

module.exports = { createNullScene, MeshBuilder };
