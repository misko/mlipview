import { createBasicScene } from "./scene.js";
import { buildBenzene } from "./benzene.js";
import { enableAtomDragging } from "./interaction.js";

const canvas = document.getElementById("renderCanvas");
const engine = new BABYLON.Engine(canvas, true, { antialias: true });

// Scene
const scene = createBasicScene(engine, canvas);

// (Optional) red debug sphere at origin – comment out once done
const debugSphere = BABYLON.MeshBuilder.CreateSphere("debugSphere", { diameter: 0.5, segments: 24 }, scene);
const debugMat = new BABYLON.StandardMaterial("debugMat", scene);
debugMat.diffuseColor = new BABYLON.Color3(1, 0, 0);
debugMat.emissiveColor = new BABYLON.Color3(0.2, 0, 0);
debugSphere.material = debugMat;

// Molecule
const { baseC, baseH, bondUnit, atoms, refreshBonds } = buildBenzene(scene, {
  carbonRingRadius: 1.9,
  hydrogenRadius: 3.1,
  bondRadius: 0.08,
  debugAlwaysActive: true
});

console.log("C:", baseC.thinInstanceCount, "H:", baseH.thinInstanceCount, "bonds:", bondUnit.thinInstanceCount);

// Enable dragging
enableAtomDragging(scene, { atoms, refreshBonds });

// WebXR (desktop-safe; shows "Enter VR" when available)
(async () => {
  try {
    const xr = await scene.createDefaultXRExperienceAsync({
      // pointer selection + teleport are enabled by default
    });
    // If you want near grabbing later:
    // xr.baseExperience.featuresManager.enableFeature(BABYLON.WebXRFeatureName.NEAR_INTERACTION, "latest");
  } catch (e) {
    console.warn("XR init skipped:", e);
  }
})();

// Render
engine.runRenderLoop(() => scene.render());
addEventListener("resize", () => engine.resize());
