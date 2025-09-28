// public/main.js - Main application entry point
import { setupScene, setupMolecule, setupBondPicking, setupGlobalFunctions } from "./setup/app-setup.js";
import { createForceField } from './forcefield/registry.js';
import { createLJForceField } from './physics_mock.js';
import { createPhysicsManager } from './physics/physics-manager.js';
import { createEnergyHUD, createStateBar } from "./ui/hud.js";
import { showMoleculeSelector } from './ui/molecule-selector.js';
import { createEnergyPlot } from "./ui/energy-plot.js";
import { createForceControls } from "./controls/force-controls.js";
import { DEBUG, dbg } from "./util/debug.js";

// DEBUG and dbg provided by util/debug.js

// Initialize application
async function initApp() {
  console.log("[app] Starting application...");
  
  // Setup 3D scene
  const { canvas, engine, scene } = await setupScene();
  
  // Setup molecule and related systems
  const { mol, atoms, state, torsion, mlip: initialMlip, forceVis } = await setupMolecule(scene);
  let mlip = initialMlip;
  
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
  
  dbg("[DEBUG] HUD elements found:", {
    energyVal: !!hud.energyVal,
    btnForces: !!hud.btnForces,
    btnLenMode: !!hud.btnLenMode,
    btnPlot: !!hud.btnPlot
  });
  
  // Create physics manager
  const physics = createPhysicsManager({ state, getMlip: () => mlip, energyPlot });
  
  // Setup UI event handlers
  setupEventHandlers(hud, stateBar, energyPlot, forceControls, state, () => mlip, physics);
  // Molecule selector button
  if (hud.btnMolecules) {
    hud.btnMolecules.onclick = () => showMoleculeSelector();
  }
  
  // Add initial data point to plot (supports async compute)
  if (energyPlot.isEnabled()) {
    const r0 = mlip.compute();
    if (r0 && typeof r0.then === 'function') {
      r0.then(res => {
        if (res && Number.isFinite(res.energy)) {
          energyPlot.addInitialDataPoint(res.energy);
        }
      }).catch(e => console.warn('[forcefield] initial compute failed', e));
    } else if (r0 && Number.isFinite(r0.energy)) {
  energyPlot.addInitialDataPoint(r0.energy);
    }
  }
  
  // Start render loop
  startRenderLoop(engine, scene, hud.energyVal, forceControls, energyPlot, state, physics);
  // Expose explicit step recorder (bond rotations will call this)
  window.recordEnergyStep = () => {
    const e = physics.energy;
    if (Number.isFinite(e) && energyPlot.isEnabled()) {
      energyPlot.recordStep(e);
    }
  };

  // Forcefield switch handlers
  function applyActiveFFStyles(kind) {
    if (hud.btnFFFair) hud.btnFFFair.style.opacity = kind === 'fairchem' ? '1' : '0.4';
    if (hud.btnFFLJ) hud.btnFFLJ.style.opacity = kind === 'lj' ? '1' : '0.4';
  }
  applyActiveFFStyles(mlip.kind || 'fairchem');
  function switchBackend(kind) {
    try {
      const created = createForceField(kind, { molecule: mol });
      mlip = created;
      console.log('[forcefield] Switched backend to', kind);
    } catch (e) {
      console.warn('[forcefield] Failed to switch to', kind, 'falling back to LJ.', e.message);
      mlip = createLJForceField(mol);
    }
    applyActiveFFStyles(mlip.kind);
    // Reset caches & plot
  physics.resetCache();
  if (energyPlot.isEnabled()) energyPlot.clear();
  }
  if (hud.btnFFFair) hud.btnFFFair.onclick = () => switchBackend('fairchem');
  if (hud.btnFFLJ) hud.btnFFLJ.onclick = () => switchBackend('lj');
  
  // Handle window resize
  addEventListener("resize", () => engine.resize());
  
  console.log("[app] Application initialized successfully");
}

function setupEventHandlers(hud, stateBar, energyPlot, forceControls, state, getMlip, physics) {
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
  dbg("[DEBUG] Plot button clicked, current state:", energyPlot.isEnabled());
    const isEnabled = energyPlot.toggle();
    hud.btnPlot.textContent = `Plot: ${isEnabled ? "ON" : "OFF"}`;
  dbg("[DEBUG] Plot toggled to:", isEnabled);
    
    if (isEnabled) {
      const r = getMlip().compute();
      if (r && typeof r.then === 'function') {
  r.then(res => res && Number.isFinite(res.energy) && energyPlot.addInitialDataPoint(res.energy));
      } else if (r && Number.isFinite(r.energy)) {
  energyPlot.addInitialDataPoint(r.energy);
      }
    }
  };
  
  // Relaxation controls
  if (hud.btnRelax) {
    hud.btnRelax.onclick = () => {
      if (physics.isRelaxationRunning()) {
        physics.stopRelaxation();
        hud.btnRelax.textContent = "Relax";
        hud.btnRelax.style.opacity = "1";
      } else {
        const options = {
          optimizer: hud.relaxOptimizer ? hud.relaxOptimizer.value : 'steepest_descent',
          stepSize: hud.relaxStepSize ? parseFloat(hud.relaxStepSize.value) : 0.01,
          maxSteps: hud.relaxMaxSteps ? parseInt(hud.relaxMaxSteps.value) : 1000,
          stepDelay: 0, // No delay - run as fast as server allows
          energyTolerance: 1e-4 // Stop when energy changes less than 1e-4
        };
        physics.startRelaxation(options);
        hud.btnRelax.textContent = "Stop";
        hud.btnRelax.style.opacity = "0.7";
      }
    };
  }

  // Toggle relaxation controls visibility with right-click
  if (hud.btnRelax && hud.relaxControls) {
    hud.btnRelax.oncontextmenu = (e) => {
      e.preventDefault();
      const isVisible = hud.relaxControls.style.display !== 'none';
      hud.relaxControls.style.display = isVisible ? 'none' : 'block';
    };
  }

  // Step size slider
  if (hud.relaxStepSize && hud.relaxStepSizeVal) {
    hud.relaxStepSize.oninput = () => {
      hud.relaxStepSizeVal.textContent = hud.relaxStepSize.value;
    };
  }
  
  // Rebuild / Export controls removed
  
  // Global clear function
  window.clearEnergyPlot = () => energyPlot.clear();
}

function updateRelaxationStatus(physics) {
  const state = physics.getRelaxationState();
  const statusEl = document.querySelector('#relaxStatus');
  const btnRelax = document.querySelector('#btnRelax');
  
  if (statusEl) {
    let statusText = `Status: ${state.status}`;
    if (state.step > 0) {
      statusText += ` | Step: ${state.step}`;
    }
    if (state.maxForce !== null) {
      statusText += ` | Max Force: ${state.maxForce.toFixed(4)}`;
    }
    if (state.energy !== null) {
      statusText += ` | Energy: ${state.energy.toFixed(6)}`;
    }
    statusEl.textContent = statusText;
  }
  
  // Update button text when relaxation finishes
  if (btnRelax && !physics.isRelaxationRunning() && state.status !== 'running') {
    btnRelax.textContent = "Relax";
    btnRelax.style.opacity = "1";
  }
}

function startRenderLoop(engine, scene, energyVal, forceControls, energyPlot, state, physics) {
  let frameCount = 0;
  let lastEnergyUpdate = 0;
  let lastForceUpdate = 0;
  // Plot updates now event-driven; remove lastPlotUpdate
  const ENERGY_UPDATE_INTERVAL = 2;
  const FORCE_UPDATE_INTERVAL = 3;
  // PLOT_UPDATE_INTERVAL removed (event-driven)
  engine.runRenderLoop(() => {
    frameCount++;
    physics.tickFrame();
    const physicsUpdated = physics.updatePhysicsCache();
    const e = physics.energy;
    const needsEnergyUpdate = physicsUpdated || (frameCount - lastEnergyUpdate >= ENERGY_UPDATE_INTERVAL);
    const needsForceUpdate = physicsUpdated || (frameCount - lastForceUpdate >= FORCE_UPDATE_INTERVAL);
  // Plot no longer updated on a timer
    if (needsEnergyUpdate) {
      if (Number.isFinite(e)) {
        energyVal.textContent = e.toFixed(3);
        lastEnergyUpdate = frameCount;
      } else if (!energyVal.textContent) {
        energyVal.textContent = 'n/a';
      }
      
      // Update relaxation status
      updateRelaxationStatus(physics);
    }
    if (needsForceUpdate && physics.forces) { forceControls.updateForces(physics.forces); lastForceUpdate = frameCount; }
    // Plot updates only happen when window.recordEnergyStep is invoked (after rotations)
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
