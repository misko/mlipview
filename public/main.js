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
  engine.runRenderLoop(() => {
    const { energy, forces } = mlip.compute();
    
    // Update energy display
    energyVal.textContent = energy.toFixed(3);
    
    // Update forces if enabled
    forceControls.updateForces(forces);
    
    // Update plot if enabled
    if (energyPlot.isEnabled()) {
      energyPlot.updateChart(energy, state);
    }
    
    scene.render();
  });
}

// Start the application
initApp().catch(console.error);
