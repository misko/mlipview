// vr/main-vr.js - VR-enabled version of the molecular viewer
import { setupVR } from './vr-setup.js';
import { createVRUI, createVRUIOnCamera } from './vr-ui.js';
// Plot disabled in VR Lite mode

// Import existing modules
import { setupScene, setupMolecule } from '../setup/app-setup.js';

export async function initVRApp() {
  console.log("[VR] 🚀 Starting VR molecular viewer...");
  console.log("[VR] Desktop app should NOT be running - checking...");
  
  // Ensure we're in VR mode, not desktop mode
  if (window.location.pathname.includes('vr.html')) {
    console.log("[VR] ✅ Confirmed we're in VR mode");
  } else {
    console.warn("[VR] ⚠️ Not in VR mode - this might cause conflicts");
  }
  
  const liteMode = true; // Force VR Lite mode only
  console.log("[VR] LITE MODE ENABLED - Performance optimizations active");
  
  try {
    console.log("[VR] Step 1: Setting up scene...");
    // Create Babylon.js app using existing setup with VR optimizations
    const { canvas, engine, scene } = await setupScene(true); // true = VR mode
    console.log("[VR] Scene setup complete:", {
      canvas: !!canvas,
      engine: !!engine,
      scene: !!scene,
      meshCount: scene.meshes.length,
      lightCount: scene.lights.length
    });
    
    console.log("[VR] Step 2: Setting up molecule...");
    // Setup molecule and related systems using existing setup
    const { mol, atoms, state, torsion, mlip, forceVis } = await setupMolecule(scene);
    console.log("[VR] Molecule setup complete:", {
      mol: !!mol,
      atoms: atoms?.length || 0,
      state: !!state,
      torsion: !!torsion,
      mlip: !!mlip,
      forceVis: !!forceVis
    });
    
    // Debug scene state before VR initialization
    console.log("[VR] Pre-VR scene debug:");
    console.log("- Scene meshes:", scene.meshes.length);
    console.log("- Scene lights:", scene.lights.length);
    console.log("- Active camera:", scene.activeCamera?.constructor.name);
    
    // Check for molecules/atoms
    const atomMeshes = scene.meshes.filter(m => m.name && m.name.includes('atom'));
    const bondMeshes = scene.meshes.filter(m => m.name && m.name.includes('bond'));
    console.log(`- Atom meshes: ${atomMeshes.length}`);
    console.log(`- Bond meshes: ${bondMeshes.length}`);
    
    console.log("[VR] Step 3: Initializing VR...");
    // Initialize VR
    const vrHelper = await setupVR(engine, scene);
    
    if (!vrHelper) {
      console.error("[VR] ❌ VR not available. setupVR returned null.");
      return null;
    }
    
    console.log("[VR] ✅ VR helper created successfully");
    
    // CRITICAL: Make the VR helper globally accessible for browser native VR button
    window.vrHelper = vrHelper;
    window.vrScene = scene;
    
  // CRITICAL: Make state, molecule, and torsion system globally available for VR interaction
    window.appState = state;
  window.vrMol = mol;
    window.vrTorsion = torsion;
    window.vrMlip = mlip;
    console.log("[VR] 🌐 VR helper and molecular systems stored globally for interaction");
    
  // Instructions panel removed (prefer 2D HUD overlay only)
    
    // Native VR button trigger is handled centrally in vr-setup.js; duplicate logic removed here
    
  console.log("[VR] Step 4: Creating VR UI...");
  // Create fullscreen HUD for non-XR; XR HUD will be created when session starts
  let vrUI = createVRUI(scene);
  let xrHud = null;
  console.log("[VR] VR UI created (energy panel)");

  // Bond rotation controls via 2D HUD overlay (bottom bar)
  try {
  const uiState = vrUI?.bond?.state || { side: 'j', step: 5, recompute: false };
  window.vrBondUI = uiState;
  let bondLabel = vrUI?.bond?.label;
  let btnMinus = vrUI?.bond?.btnMinus;
  let btnPlus  = vrUI?.bond?.btnPlus;

  // Wire clicks to rotation (with debug)
  btnMinus?.onPointerUpObservable.add(() => { console.log('[HUD] (-) pressed'); rotateSelectedBond(scene, -1); });
  btnPlus?.onPointerUpObservable.add(() => { console.log('[HUD] (+) pressed'); rotateSelectedBond(scene, +1); });

    // Update label based on selection (bar stays visible to make HUD obvious)
    const __updateBondHud = () => {
      const sel = window.vrSelection;
      const label = bondLabel || xrHud?.bond?.label;
      if (!label) return;
      if (sel && sel.kind === 'bond') label.text = `Bond ${sel.data.i}-${sel.data.j}`; else label.text = 'Select a bond to rotate';
    };
    scene.onBeforeRenderObservable.add(__updateBondHud);
    // Build a camera-anchored HUD for XR session
    vrHelper.baseExperience.sessionManager.onXRSessionInit.add(() => {
      try {
        const xrCam = vrHelper?.baseExperience?.camera || scene.activeCamera;
        if (!xrCam) return;
        if (xrHud && xrHud.advancedTexture) { try { xrHud.advancedTexture.dispose(); } catch {} }
        if (xrHud && xrHud.rootMesh) { try { xrHud.rootMesh.dispose(); } catch {} }
        xrHud = createVRUIOnCamera(scene, xrCam);
        // Rebind button handlers to XR HUD controls
        btnMinus = xrHud?.bond?.btnMinus || btnMinus;
        btnPlus  = xrHud?.bond?.btnPlus  || btnPlus;
        bondLabel = xrHud?.bond?.label   || bondLabel;
        // Use XR HUD state for torsion controls going forward
        try {
          const prev = window.vrBondUI || { side: 'j', step: 5, recompute: false };
          if (xrHud?.bond?.state) {
            xrHud.bond.state.side = prev.side || 'j';
            xrHud.bond.state.step = Math.abs(prev.step || 5);
            xrHud.bond.state.recompute = !!prev.recompute;
            window.vrBondUI = xrHud.bond.state;
            // Sync visible labels for side/recompute
            if (xrHud.bond.btnSide?.textBlock) xrHud.bond.btnSide.textBlock.text = `Side: ${window.vrBondUI.side}`;
            if (xrHud.bond.btnRec?.textBlock) xrHud.bond.btnRec.textBlock.text = `recompute: ${window.vrBondUI.recompute ? 'on' : 'off'}`;
          }
        } catch {}
        btnMinus?.onPointerUpObservable.add(() => { console.log('[HUD XR] (-) pressed'); rotateSelectedBond(scene, -1); });
        btnPlus?.onPointerUpObservable.add(() => { console.log('[HUD XR] (+) pressed'); rotateSelectedBond(scene, +1); });
      } catch {}
      console.log('[VR] Created camera-anchored HUD for XR session');
    });
    vrHelper.baseExperience.sessionManager.onXRSessionEnded.add(() => {
      try { if (xrHud?.advancedTexture?.dispose) xrHud.advancedTexture.dispose(); } catch {}
      try { if (xrHud?.rootMesh?.dispose) xrHud.rootMesh.dispose(); } catch {}
      xrHud = null;
      // fullscreen HUD remains available outside XR
    });
  } catch (e) {
    console.warn('[VR] Overlay bond controls failed:', e?.message);
  }
    
    // VR-specific state
    // Lite mode: remove advanced UI control wiring (forces, plot, export)
    let vrUIVisible = true;
    
    // Update energy display in VR
    function updateVREnergy() {
      const { energy } = mlip.compute();
      vrUI.energyValue.text = energy.toFixed(3);
      
      // Plot disabled in lite mode
    }
    
    // Render loop with cached physics results
    function startVRRenderLoop() {
      let frameCount = 0;
      let lastEnergyUpdate = 0;
  // Removed: force and plot update trackers
      
      // Cache physics results - only recompute when molecular structure changes
      let cachedEnergy = null;
  // Removed: cached forces (not used in lite mode)
      let lastChangeCounter = -1;
      
      // Function to check if molecular structure has changed
      function hasStructureChanged() {
        if (mol && mol.changeCounter !== lastChangeCounter) {
          lastChangeCounter = mol.changeCounter;
          return true;
        }
        return false;
      }
      
      // Function to update cached physics if needed
      function updatePhysicsCache() {
        if (hasStructureChanged() || cachedEnergy === null) {
          console.log("[VR] Molecular structure changed - recomputing physics");
          const result = mlip.compute();
          cachedEnergy = result.energy;
          return true;
        }
        return false;
      }
      
      // Throttle UI updates (less aggressive now since physics is cached)
  const ENERGY_UPDATE_INTERVAL = 5;   // Update energy display every 5 frames in lite mode
      
      engine.runRenderLoop(() => {
        try {
          frameCount++;
          
          // Check if we need to recompute physics (only when structure changes)
          const physicsUpdated = updatePhysicsCache();
          
          // Determine what UI elements need updating
          const needsEnergyUpdate = physicsUpdated || (frameCount - lastEnergyUpdate >= ENERGY_UPDATE_INTERVAL);
          // Lite mode: only energy update used
          
          // Update VR energy display (using cached energy)
          if (needsEnergyUpdate && vrUI && vrUI.energyValue && cachedEnergy !== null) {
            vrUI.energyValue.text = cachedEnergy.toFixed(3);
            lastEnergyUpdate = frameCount;
          }

          // Handle LEFT-stick bond rotation (up/down) when a bond is selected
          try {
            const sel = window.vrSelection;
            if (sel && sel.kind === 'bond') {
              const helper = window.vrHelper;
              const controllers = helper?.input?.controllers || [];
              const deadzone = 0.25; // require intent
              // Throttle to avoid too-fast repeats
              if (!window.__vrBondStickNext) window.__vrBondStickNext = 0;
              const now = performance.now();
              const repeatMs = 130; // step roughly ~7.5/s when held
              let leftY = 0;
              // Prefer LEFT-hand motionController thumbstick component if present
              let leftCtrl = controllers.find(c => c?.inputSource?.handedness === 'left');
              const thumb = leftCtrl?.motionController?.getComponent?.('xr-standard-thumbstick');
              if (thumb && thumb.axes && typeof thumb.axes.y === 'number') {
                leftY = thumb.axes.y;
                if (window && window.vrLogSticks && typeof console !== 'undefined' && console.debug) {
                  console.debug('[VR Stick] Left thumbstick axes (motionController):', { x: thumb.axes.x, y: thumb.axes.y });
                }
              } else if (leftCtrl?.inputSource?.gamepad?.axes) {
                const gp = leftCtrl.inputSource.gamepad;
                leftY = gp.axes.length >= 2 ? gp.axes[1] : 0;
                if (window && window.vrLogSticks && typeof console !== 'undefined' && console.debug) {
                  console.debug('[VR Stick] Left controller axes (gamepad):', JSON.stringify(Array.from(gp.axes)));
                }
              } else {
                // Fallback: pick any controller with strongest vertical axis
                for (const c of controllers) {
                  const gp = c.inputSource?.gamepad;
                  if (!gp || !gp.axes) continue;
                  const ay = gp.axes.length >= 2 ? gp.axes[1] : 0;
                  if (Math.abs(ay) > Math.abs(leftY)) leftY = ay;
                }
                if (window && window.vrLogSticks && typeof console !== 'undefined' && console.debug) {
                  const axDump = controllers.map(c => Array.from(c.inputSource?.gamepad?.axes || []) );
                  console.debug('[VR Stick] Fallback axes dump:', JSON.stringify(axDump));
                }
              }
              if (Math.abs(leftY) > deadzone && now >= window.__vrBondStickNext) {
                const sign = leftY < 0 ? +1 : -1; // up is negative Y: rotate +, down is positive Y: rotate -
                if (window && window.vrLogSticks && typeof console !== 'undefined') {
                  console.log('[VR Stick] Rotate via left stick', { leftY: +leftY.toFixed(3), sign });
                }
                rotateSelectedBond(scene, sign);
                window.__vrBondStickNext = now + repeatMs;
              } else if (Math.abs(leftY) <= deadzone) {
                // Reset throttle quickly when stick returns to center
                window.__vrBondStickNext = 0;
              }
            }
          } catch {}
          
          // Update forces if enabled (using cached forces) - skip in lite mode
          // Forces disabled in lite mode
          
          // Update plot if enabled (using cached energy) - skip in lite mode
          // Plot disabled in lite mode
          
          scene.render();
        } catch (error) {
          console.error('[VR Render Loop] Error:', error);
          // Don't break the render loop, just log the error
        }
      });
    }
    
    // Initialize default state
  // Lite mode: forces/plot disabled
  updateVREnergy();
    
    startVRRenderLoop();
    
    // Handle window resize
    addEventListener("resize", () => engine.resize());
    
    console.log("[VR] VR molecular viewer initialized successfully!");
    
    // Final debug check
    console.log("[VR] Final scene state:");
    console.log(`- Total meshes: ${scene.meshes.length}`);
    console.log(`- Visible meshes: ${scene.meshes.filter(m => m.isVisible).length}`);
    console.log(`- Total lights: ${scene.lights.length}`);
    console.log(`- Enabled lights: ${scene.lights.filter(l => l.isEnabled()).length}`);
    console.log(`- Camera position: ${scene.activeCamera?.position}`);
    
    // Expose debug info globally for Quest console access
    window.vrDebug = {
      scene,
      engine,
      vrHelper,
      mol,
      state,
  forceVis,
  vrUI,
      mlip, // Add MLIP for debugging
      debugInfo: () => {
        console.log("=== VR DEBUG INFO ===");
        console.log("Scene meshes:", scene.meshes.map(m => `${m.name} (visible: ${m.isVisible})`));
        console.log("Scene lights:", scene.lights.map(l => `${l.constructor.name} (intensity: ${l.intensity}, enabled: ${l.isEnabled()})`));
        console.log("Camera:", scene.activeCamera);
        console.log("VR Session active:", vrHelper.baseExperience.sessionManager.currentSession !== null);
        console.log("Molecule change counter:", mol.changeCounter);
        console.log("=== END DEBUG ===");
      }
    };
    
    return {
      engine,
      scene,
      mol,
      state,
      vrHelper,
      vrUI
    };
    
  } catch (error) {
    console.error("[VR] Failed to initialize VR app:", error);
    return null;
  }
}

// Rotate the currently selected bond by a fixed step via torsion controller
function rotateSelectedBond(scene, sign) {
  const sel = window.vrSelection;
  const torsion = window.vrTorsion;
  if (!sel) { console.warn('[VR] rotateSelectedBond called but no selection'); return; }
  if (sel.kind !== 'bond') { console.warn('[VR] rotateSelectedBond called but selection is not a bond:', sel.kind); return; }
  if (!torsion) { console.warn('[VR] rotateSelectedBond called but torsion controller missing'); return; }
  const ui = window.vrBondUI || { side: 'j', step: 5, recompute: false };
  const step = Math.abs(ui.step || 5);
  try {
    if (typeof console !== 'undefined') {
      console.log('[VR] rotateSelectedBond', { i: sel.data.i, j: sel.data.j, side: ui.side, step, sign, recompute: !!ui.recompute });
    }
    torsion.rotateAroundBond({ i: sel.data.i, j: sel.data.j, side: ui.side || 'j', angleDeg: sign * step, recompute: !!ui.recompute });
    if (typeof window.vrMol?.markChanged === 'function') window.vrMol.markChanged();
    // One-shot recompute like desktop: reset flag after use
    if (ui.recompute) ui.recompute = false;
    // Keep bond selected: do not clear window.vrSelection here.
    // Re-emit label update by touching a flag (safe no-op)
    try { scene._lastHudTick = (scene._lastHudTick || 0) + 1; } catch {}
  } catch (e) {
    console.warn('[VR] rotateSelectedBond failed:', e?.message);
  }
}

function createVRInstructions(scene) {
  try {
    // Create a 3D text panel with VR instructions
    const manager = new BABYLON.GUI.GUI3DManager(scene);
    
    // Create a panel
    const panel = new BABYLON.GUI.CylinderPanel();
    panel.margin = 0.2;
    
    // Position the panel in front of the user
    panel.position.z = 2;
    panel.position.y = 1.5;
    panel.radius = 3;
    
    manager.addControl(panel);
    
    const instructions = [
      "🎮 VR Controls (Lite):",
      "",
      "✋ Grab box: Move/rotate the molecule",
      "🤲 Two hands: Pinch to scale",
      " Tap trigger: Select bond/atom",
      " ⬆️⬇️ Left stick: Rotate selected bond",
      "👉 Hold trigger + rotate wrist (fallback)",
      "👀 Move your head to look around",
      "",
      "🧪 Molecule: ROY Crystal Structure"
    ];
    
    instructions.forEach((text, index) => {
      const button = new BABYLON.GUI.HolographicButton("instruction_" + index);
      panel.addControl(button);
      
      button.text = text;
      button.fontSize = 24;
      button.color = "white";
      button.background = "rgba(0, 100, 200, 0.8)";
      
      // Make it non-interactive for pure display
      button.isPointerBlocker = false;
      button.isEnabled = false;
    });
    
    // Make the panel fade out after 10 seconds
    setTimeout(() => {
      let alpha = 1;
      const fadeInterval = setInterval(() => {
        alpha -= 0.05;
        panel.scaling = new BABYLON.Vector3(alpha, alpha, alpha);
        
        if (alpha <= 0) {
          clearInterval(fadeInterval);
          manager.removeControl(panel);
          panel.dispose();
        }
      }, 100);
    }, 10000);
    
    console.log("[VR] Instructions panel created and will auto-hide in 10 seconds");
    
  } catch (error) {
    console.warn("[VR] Could not create instructions panel:", error.message);
  }
}

// showVRTextPanel removed (unused)
