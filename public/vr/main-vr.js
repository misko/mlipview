// vr/main-vr.js - VR-enabled version of the molecular viewer
import { setupVR } from './vr-setup.js';
import { createVRUI } from './vr-ui.js';
import { create3DEnergyPlot } from './vr-plot.js';

// Import existing modules
import { setupScene, setupMolecule } from '../setup/app-setup.js';

export async function initVRApp() {
  console.log("[VR] Starting VR molecular viewer...");
  
  try {
    // Create Babylon.js app using existing setup
    const { canvas, engine, scene } = await setupScene();
    
    // Setup molecule and related systems using existing setup
    const { mol, atoms, state, torsion, mlip, forceVis } = await setupMolecule(scene);
    
    // Initialize VR
    const vrHelper = await setupVR(engine, scene);
    
    if (!vrHelper) {
      console.error("VR not available. Falling back to desktop mode.");
      return null;
    }
    
    // Create VR UI components
    const vrUI = createVRUI(scene);
    const vr3DPlot = create3DEnergyPlot(scene);
    
    // VR-specific state
    let vrUIVisible = true;
    let vrPlotVisible = true;
    let vrForcesVisible = true;
    
    // Connect VR UI buttons to functionality
    vrUI.buttons.forces.onPointerClickObservable.add(() => {
      vrForcesVisible = !vrForcesVisible;
      vrUI.buttons.forces.textBlock.text = `Forces: ${vrForcesVisible ? 'ON' : 'OFF'}`;
      forceVis.setEnabled(vrForcesVisible);
    });
    
    vrUI.buttons.plot.onPointerClickObservable.add(() => {
      vrPlotVisible = !vrPlotVisible;
      vrUI.buttons.plot.textBlock.text = `Plot: ${vrPlotVisible ? 'ON' : 'OFF'}`;
      vr3DPlot.toggle(vrPlotVisible);
    });
    
    vrUI.buttons.length.onPointerClickObservable.add(() => {
      // Toggle between normalized and true force lengths
      const currentMode = forceVis.getLengthMode();
      const newMode = currentMode === "normalized" ? "true" : "normalized";
      forceVis.setLengthMode(newMode);
      vrUI.buttons.length.textBlock.text = `Length: ${newMode === 'true' ? 'True' : 'Normalized'}`;
    });
    
    vrUI.buttons.rebuild.onPointerClickObservable.add(() => {
      console.log("[VR] Rebuilding from rotations...");
      state.recomputeAndCommit();
      updateVREnergy();
    });
    
    vrUI.buttons.export.onPointerClickObservable.add(() => {
      const xyz = state.exportXYZ("VR_EXPORT");
      console.log("[VR] XYZ Export:", xyz);
      showVRTextPanel(xyz);
    });
    
    // Update energy display in VR
    function updateVREnergy() {
      const { energy } = mlip.compute();
      vrUI.energyValue.text = energy.toFixed(3);
      
      if (vrPlotVisible) {
        vr3DPlot.addDataPoint(state.rotations.length, energy);
      }
    }
    
    // Render loop
    function startVRRenderLoop() {
      engine.runRenderLoop(() => {
        const { energy, forces } = mlip.compute();
        
        // Update VR energy display
        vrUI.energyValue.text = energy.toFixed(3);
        
        // Update forces if enabled
        if (vrForcesVisible) {
          forceVis.updateForces(forces);
        }
        
        // Update plot if enabled
        if (vrPlotVisible) {
          vr3DPlot.addDataPoint(state.rotations.length, energy);
        }
        
        scene.render();
      });
    }
    
    // Initialize default state
    forceVis.setEnabled(vrForcesVisible);
    vr3DPlot.toggle(vrPlotVisible);
    updateVREnergy();
    
    startVRRenderLoop();
    
    // Handle window resize
    addEventListener("resize", () => engine.resize());
    
    console.log("[VR] VR molecular viewer initialized successfully!");
    
    return {
      engine,
      scene,
      mol,
      state,
      vrHelper,
      vrUI,
      vr3DPlot
    };
    
  } catch (error) {
    console.error("[VR] Failed to initialize VR app:", error);
    return null;
  }
}

function showVRTextPanel(text) {
  console.log("VR Text Panel would show:", text);
}
