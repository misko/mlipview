// setup/app-setup.js - Main application setup
import { createBasicScene } from "../scene.js";
import { enableAtomDragging } from "../interaction.js";
import { loadBondsSpec } from "../loader_bonds.js";
import { createTorsionController } from "../torsion.js";
import { createMockMLIP } from "../physics_mock.js";
import { createForceRenderer } from "../forces.js";
import { enableBondPicking } from "../bond_pick.js";
import { createStateStore } from "../state.js";
import { buildDefault } from "../molecules/molecule-loader.js";

export async function setupScene(vrMode = false) {
  const canvas = document.getElementById("renderCanvas");
  const engine = new BABYLON.Engine(canvas, true, { 
    antialias: !vrMode, // Disable antialiasing in VR for performance
    stencil: !vrMode,   // Disable stencil buffer in VR if not needed
    preserveDrawingBuffer: false,
    powerPreference: "high-performance"
  });
  
  const scene = createBasicScene(engine, canvas);
  
  // VR performance optimizations
  if (vrMode) {
    console.log("[VR] Applying VR performance optimizations...");
    
    // Reduce quality settings for VR
    scene.skipPointerMovePicking = true;
    scene.skipFrustumClipping = false; // Keep frustum culling for performance
    scene.autoClear = true;
    scene.autoClearDepthAndStencil = true;
    
    // Optimize rendering pipeline
    scene.getEngine().setHardwareScalingLevel(1.0);
    
    // ENHANCED LIGHTING FOR VR - Don't reduce lighting intensity, instead increase it
    scene.lights.forEach(light => {
      if (light instanceof BABYLON.HemisphericLight) {
        light.intensity = Math.max(light.intensity, 1.5); // Ensure bright hemisphere light
        console.log("[VR] Enhanced hemisphere light intensity to", light.intensity);
      } else if (light instanceof BABYLON.DirectionalLight) {
        light.intensity = Math.max(light.intensity, 1.0); // Ensure bright directional light
        console.log("[VR] Enhanced directional light intensity to", light.intensity);
      }
    });
    
    // Add additional VR lighting to prevent black screen
    const vrBoostLight = new BABYLON.HemisphericLight("vrBoost", new BABYLON.Vector3(0, 1, 0), scene);
    vrBoostLight.intensity = 1.2;
    vrBoostLight.diffuse = new BABYLON.Color3(1, 1, 1);
    console.log("[VR] Added VR boost lighting");
    
    console.log("[VR] Performance optimizations applied with enhanced lighting");
  }
  
  return { canvas, engine, scene };
}

export async function setupMolecule(scene) {
  const mol = await buildDefault(scene);
  const { atoms, refreshBonds } = mol;
  
  // Add change tracking to the molecule
  mol.changeCounter = 0;
  mol.markChanged = function() {
    this.changeCounter++;
    console.log(`[Molecule] Structure changed (change #${this.changeCounter})`);
  };
  
  // Persistent state (records torsions & can reconstruct/export)
  const state = createStateStore(mol);
  
  // Drag interaction
  enableAtomDragging(scene, { atoms, refreshBonds, molecule: mol });
  
  // Torsion controller (records into state)
  const torsion = createTorsionController(mol, state);
  
  // Mock MLIP + force visualization
  const mlip = createMockMLIP(mol);
  const forceVis = createForceRenderer(scene, atoms, {
    color: new BABYLON.Color3(1, 0.2, 0.9)
  });
  
  return { mol, atoms, state, torsion, mlip, forceVis };
}

export async function setupBondPicking(scene, mol, torsion) {
  try {
    const rotatable = await loadBondsSpec("./molecules/roy.bonds");
    console.log("[torsion] Loaded rotatable spec:", rotatable);
    
    enableBondPicking(scene, {
      molecule: mol,
      torsion,
      rotatableSpec: rotatable,
      strict: false
    });
    
    return rotatable;
  } catch (e) {
    console.warn("[torsion] No roy.bonds found (optional).", e);
    
    enableBondPicking(scene, {
      molecule: mol,
      torsion,
      rotatableSpec: [],
      strict: false
    });
    
    return [];
  }
}

export function setupGlobalFunctions(state, torsion) {
  // Global functions for console access
  window.stateJson = () => state.getStateJSON();
  window.stateExportXYZ = (name) => state.exportXYZ(name);
  window.rotateBond = (i, j, side = "j", angleDeg = 5, recompute = false) =>
    torsion.rotateAroundBond({ i, j, side, angleDeg, recompute });
}
