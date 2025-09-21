// public/main.js - Main application entry point
import { setupScene, setupMolecule, setupBondPicking, setupGlobalFunctions } from "./setup/app-setup.js";
import { createEnergyHUD, createStateBar } from "./ui/hud.js";
import { createEnergyPlot } from "./ui/energy-plot.js";
import { createForceControls } from "./controls/force-controls.js";

// Initialize application
async function initApp() {
  console.log("[app] Starting application...");
  
  // Setup 3D scene
  const { canvas, engine, scene } = await setupScene();
  
  // Setup molecule and related systems
  const { mol, atoms, state, torsion, mlip, forceVis } = await setupMolecule(scene);
  
  // Setup bond picking
  await setupBondPicking(scene, mol, torsion);
  
  // Setup global functions
  setupGlobalFunctions(state, torsion);
  
  // Make state globally available for plot clearing
  window.appState = state;
  
  // Create UI components
  const hud = createEnergyHUD();
  const stateBar = createStateBar();
  const energyPlot = createEnergyPlot();
  const forceControls = createForceControls(forceVis);
  
  console.log("[DEBUG] HUD elements found:", {
    energyVal: !!hud.energyVal,
    btnForces: !!hud.btnForces,
    btnLenMode: !!hud.btnLenMode,
    btnPlot: !!hud.btnPlot
  });
  
  // Setup UI event handlers
  setupEventHandlers(hud, stateBar, energyPlot, forceControls, state, mlip);
  
  // Add initial data point to plot
  if (energyPlot.isEnabled()) {
    const { energy } = mlip.compute();
    energyPlot.addInitialDataPoint(energy, state);
  }
  
  // Start render loop
  startRenderLoop(engine, scene, mlip, hud.energyVal, forceControls, energyPlot, state);
  
  // Handle window resize
  addEventListener("resize", () => engine.resize());
  
  console.log("[app] Application initialized successfully");
}

function setupEventHandlers(hud, stateBar, energyPlot, forceControls, state, mlip) {
  // Forces button
  hud.btnForces.onclick = () => {
    const isEnabled = forceControls.toggleForces();
    hud.btnForces.textContent = `Forces: ${isEnabled ? "ON" : "OFF"}`;
  };
  
  // Length mode button
  hud.btnLenMode.onclick = () => {
    const mode = forceControls.toggleLengthMode();
    hud.btnLenMode.textContent = `Length: ${mode === "normalized" ? "Normalized" : "True"}`;
  };
  
  // Plot button
  hud.btnPlot.onclick = () => {
    console.log("[DEBUG] Plot button clicked, current state:", energyPlot.isEnabled());
    const isEnabled = energyPlot.toggle();
    hud.btnPlot.textContent = `Plot: ${isEnabled ? "ON" : "OFF"}`;
    console.log("[DEBUG] Plot toggled to:", isEnabled);
    
    if (isEnabled) {
      const { energy } = mlip.compute();
      energyPlot.addInitialDataPoint(energy, state);
    }
  };
  
  // State bar buttons
  stateBar.btnRebuild.onclick = () => {
    state.recomputeAndCommit();
  };
  
  stateBar.btnExport.onclick = () => {
    const xyz = state.exportXYZ("ROY_from_rotations");
    const blob = new Blob([xyz], { type: "text/plain" });
    const url = URL.createObjectURL(blob);
    const a = document.createElement("a");
    a.href = url;
    a.download = "reconstructed.xyz";
    a.click();
    URL.revokeObjectURL(url);
  };
  
  // Global clear function
  window.clearEnergyPlot = () => energyPlot.clear(state);
}

function startRenderLoop(engine, scene, mlip, energyVal, forceControls, energyPlot, state) {
  // Cache physics results - only recompute when molecular structure changes
  let cachedEnergy = null;
  let cachedForces = null;
  let lastChangeCounter = -1;
  let frameCount = 0;
  let lastEnergyUpdate = 0;
  let lastForceUpdate = 0;
  let lastPlotUpdate = 0;
  
  // Function to check if molecular structure has changed
  function hasStructureChanged() {
    const molecule = state.molecule;
    if (molecule && molecule.changeCounter !== lastChangeCounter) {
      lastChangeCounter = molecule.changeCounter;
      return true;
    }
    return false;
  }
  
  // Function to update cached physics if needed
  function updatePhysicsCache() {
    if (hasStructureChanged() || cachedEnergy === null) {
      console.log("[Desktop] Molecular structure changed - recomputing physics");
      const result = mlip.compute();
      cachedEnergy = result.energy;
      cachedForces = result.forces;
      return true;
    }
    return false;
  }
  
  // Throttle UI updates
  const ENERGY_UPDATE_INTERVAL = 2;  // Update energy display every 2 frames
  const FORCE_UPDATE_INTERVAL = 3;   // Update forces every 3 frames  
  const PLOT_UPDATE_INTERVAL = 5;    // Update plot every 5 frames
  
  engine.runRenderLoop(() => {
    frameCount++;
    
    // Check if we need to recompute physics (only when structure changes)
    const physicsUpdated = updatePhysicsCache();
    
    // Determine what UI elements need updating
    const needsEnergyUpdate = physicsUpdated || (frameCount - lastEnergyUpdate >= ENERGY_UPDATE_INTERVAL);
    const needsForceUpdate = physicsUpdated || (frameCount - lastForceUpdate >= FORCE_UPDATE_INTERVAL);
    const needsPlotUpdate = physicsUpdated || (frameCount - lastPlotUpdate >= PLOT_UPDATE_INTERVAL);
    
    // Update energy display (using cached energy)
    if (needsEnergyUpdate && cachedEnergy !== null) {
      energyVal.textContent = cachedEnergy.toFixed(3);
      lastEnergyUpdate = frameCount;
    }
    
    // Update forces if enabled (using cached forces)
    if (needsForceUpdate && cachedForces) {
      forceControls.updateForces(cachedForces);
      lastForceUpdate = frameCount;
    }
    
    // Update plot if enabled (using cached energy)
    if (needsPlotUpdate && energyPlot.isEnabled() && cachedEnergy !== null) {
      energyPlot.updateChart(cachedEnergy, state);
      lastPlotUpdate = frameCount;
    }
    
    scene.render();
  });
}

// Export for use in other modules (like VR mode)
export { initApp };

// Only start the application if we're not in VR mode
// VR mode will be handled by vr.html and main-vr.js
if (!window.location.pathname.includes('vr.html')) {
  // Start the application
  initApp().catch(console.error);
}
