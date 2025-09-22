// vr/main-vr.js - VR-enabled version of the molecular viewer
import { setupVR } from './vr-setup.js';
import { createVRUI, createVRUIOnCamera } from './vr-ui.js';
// Plot disabled in VR Lite mode

// Import existing modules
import { setupScene, setupMolecule } from '../setup/app-setup.js';

// Lightweight logging guards (enable by setting window.vrLog = true in console)
function vrLog() {
  try { if (window && window.vrLog) console.log.apply(console, arguments); } catch (_) {}
}
function vrDebug() {
  try {
    if (window && window.vrLog) {
      (console.debug || console.log).apply(console, arguments);
    }
  } catch (_) {}
}

export async function initVRApp() {
  vrLog("[VR] 🚀 Starting VR molecular viewer...");
  vrLog("[VR] Desktop app should NOT be running - checking...");
  
  // Ensure we're in VR mode, not desktop mode
  if (window.location.pathname.includes('vr.html')) {
    vrLog("[VR] ✅ Confirmed we're in VR mode");
  } else {
    console.warn("[VR] ⚠️ Not in VR mode - this might cause conflicts");
  }
  
  vrLog("[VR] LITE MODE ENABLED - Performance optimizations active");
  
  try {
    vrLog("[VR] Step 1: Setting up scene...");
    // Create Babylon.js app using existing setup with VR optimizations
    const { canvas, engine, scene } = await setupScene(true); // true = VR mode
    vrDebug("[VR] Scene setup complete:", {
      canvas: !!canvas,
      engine: !!engine,
      scene: !!scene,
      meshCount: scene.meshes.length,
      lightCount: scene.lights.length
    });
    
    vrLog("[VR] Step 2: Setting up molecule...");
    // Setup molecule and related systems using existing setup
    const { mol, atoms, state, torsion, mlip, forceVis } = await setupMolecule(scene);
    vrDebug("[VR] Molecule setup complete:", {
      mol: !!mol,
      atoms: atoms?.length || 0,
      state: !!state,
      torsion: !!torsion,
      mlip: !!mlip,
      forceVis: !!forceVis
    });
    
    // Debug scene state before VR initialization
    vrDebug("[VR] Pre-VR scene debug:");
    vrDebug("- Scene meshes:", scene.meshes.length);
    vrDebug("- Scene lights:", scene.lights.length);
    vrDebug("- Active camera:", scene.activeCamera?.constructor.name);
    
    // Check for molecules/atoms
    const atomMeshes = scene.meshes.filter(m => m.name && m.name.includes('atom'));
    const bondMeshes = scene.meshes.filter(m => m.name && m.name.includes('bond'));
  vrDebug(`- Atom meshes: ${atomMeshes.length}`);
  vrDebug(`- Bond meshes: ${bondMeshes.length}`);
    
  vrLog("[VR] Step 3: Initializing VR...");
    // Initialize VR
    const vrHelper = await setupVR(engine, scene);
    
    if (!vrHelper) {
      console.error("[VR] ❌ VR not available. setupVR returned null.");
      return null;
    }
    
  vrLog("[VR] ✅ VR helper created successfully");
    
    // CRITICAL: Make the VR helper globally accessible for browser native VR button
    window.vrHelper = vrHelper;
    window.vrScene = scene;
    
  // CRITICAL: Make state, molecule, and torsion system globally available for VR interaction
    window.appState = state;
  window.vrMol = mol;
    window.vrTorsion = torsion;
    window.vrMlip = mlip;
  vrLog("[VR] 🌐 VR helper and molecular systems stored globally for interaction");
    
  // Native VR UI is handled centrally in vr-setup.js
    
  vrLog("[VR] Step 4: Creating VR UI...");
  // Create fullscreen HUD for non-XR; XR HUD will be created when session starts
  let vrUI = createVRUI(scene);
  let xrHud = null;
  vrLog("[VR] VR UI created (energy panel)");

  // Bond rotation controls via 2D HUD overlay
  try {
  const uiState = vrUI?.bond?.state || { side: 'j', step: 5, recompute: false };
  window.vrBondUI = uiState;
  // No bond label in Lite HUD
  let btnMinus = vrUI?.bond?.btnMinus;
  let btnPlus  = vrUI?.bond?.btnPlus;

  // Wire clicks to rotation (with debug)
  btnMinus?.onPointerUpObservable.add(() => { vrDebug('[HUD] (-) pressed'); rotateSelectedBond(scene, -1); try { const e = mlip.compute()?.energy; if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e); } catch {} });
  btnPlus?.onPointerUpObservable.add(() => { vrDebug('[HUD] (+) pressed'); rotateSelectedBond(scene, +1); try { const e = mlip.compute()?.energy; if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e); } catch {} });

  // Lite HUD: no bond label
    // Build a camera-anchored HUD for XR session
    vrHelper.baseExperience.sessionManager.onXRSessionInit.add(() => {
      try {
        const xrCam = vrHelper?.baseExperience?.camera || scene.activeCamera;
        if (!xrCam) return;
        // Hide 2D overlay during XR to avoid duplicate energy panel in HMD
        try { if (vrUI?.advancedTexture) vrUI.advancedTexture.rootContainer.isVisible = false; } catch {}
  // Dispose any prior XR HUD instance
  try { xrHud?.dispose && xrHud.dispose(); } catch {}
  xrHud = createVRUIOnCamera(scene, xrCam);
        // Rebind button handlers to XR HUD controls
        btnMinus = xrHud?.bond?.btnMinus || btnMinus;
        btnPlus  = xrHud?.bond?.btnPlus  || btnPlus;
        // Use XR HUD state for torsion controls going forward
        try {
          const prev = window.vrBondUI || { side: 'j', step: 5, recompute: false };
          if (xrHud?.bond?.state) {
            xrHud.bond.state.side = prev.side || 'j';
            xrHud.bond.state.step = Math.abs(prev.step || 5);
            xrHud.bond.state.recompute = !!prev.recompute;
            window.vrBondUI = xrHud.bond.state;
            // Sync visible labels for side/recompute
            if (xrHud.bond.btnSide?.textBlock) {
              const lbl = (window.vrBondUI.side === 'j') ? '(i,j)' : '(j,i)';
              xrHud.bond.btnSide.textBlock.text = lbl;
            }
            if (xrHud.bond.btnRec?.textBlock) xrHud.bond.btnRec.textBlock.text = `recompute: ${window.vrBondUI.recompute ? 'on' : 'off'}`;
          }
        } catch {}
  btnMinus?.onPointerUpObservable.add(() => {
    vrDebug('[HUD XR] (-) pressed');
    rotateSelectedBond(scene, -1);
    try {
      const e = mlip.compute()?.energy;
      if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e);
      if (xrHud?.energyValue && isFinite(e)) xrHud.energyValue.text = e.toFixed(3);
    } catch {}
  });
  btnPlus?.onPointerUpObservable.add(() => {
    vrDebug('[HUD XR] (+) pressed');
    rotateSelectedBond(scene, +1);
    try {
      const e = mlip.compute()?.energy;
      if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e);
      if (xrHud?.energyValue && isFinite(e)) xrHud.energyValue.text = e.toFixed(3);
    } catch {}
  });
  // Seed plot with current energy (if present)
        try {
          const e0 = mlip.compute()?.energy;
          if (xrHud?.plot && isFinite(e0)) xrHud.plot.addPoint(e0);
        } catch {}
      } catch {}
      vrLog('[VR] Created camera-anchored HUD for XR session');
    });
    vrHelper.baseExperience.sessionManager.onXRSessionEnded.add(() => {
      try { xrHud?.dispose && xrHud.dispose(); } catch {}
      xrHud = null;
      // Restore 2D overlay when exiting XR
      try { if (vrUI?.advancedTexture) vrUI.advancedTexture.rootContainer.isVisible = true; } catch {}
      // fullscreen HUD remains available outside XR
    });
  } catch (e) {
    console.warn('[VR] Overlay bond controls failed:', e?.message);
  }
    
    // VR-specific state: advanced features disabled in Lite mode
    
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
          vrDebug("[VR] Molecular structure changed - recomputing physics");
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
          if (needsEnergyUpdate && cachedEnergy !== null) {
            if (vrUI && vrUI.energyValue) {
              vrUI.energyValue.text = cachedEnergy.toFixed(3);
            }
            if (xrHud && xrHud.energyValue) {
              xrHud.energyValue.text = cachedEnergy.toFixed(3);
            }
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
                  vrLog('[VR Stick] Rotate via left stick', { leftY: +leftY.toFixed(3), sign });
                }
                rotateSelectedBond(scene, sign);
                // Push a point to XR HUD plot if present (energy vs step)
                try {
                  const e = (cachedEnergy !== null) ? cachedEnergy : (mlip.compute()?.energy ?? null);
                  if (xrHud?.plot && e !== null && isFinite(e)) {
                    xrHud.plot.addPoint(e);
                  }
                } catch {}
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
  // Seed energy once before render loop (Lite mode)
  try {
    const { energy } = mlip.compute();
    if (vrUI?.energyValue) vrUI.energyValue.text = energy.toFixed(3);
  } catch {}
    
  startVRRenderLoop();
    
    // Handle window resize
  addEventListener("resize", () => engine.resize());
    
  vrLog("[VR] VR molecular viewer initialized successfully!");
    
    // Final debug check
  vrDebug("[VR] Final scene state:");
  vrDebug(`- Total meshes: ${scene.meshes.length}`);
  vrDebug(`- Visible meshes: ${scene.meshes.filter(m => m.isVisible).length}`);
  vrDebug(`- Total lights: ${scene.lights.length}`);
  vrDebug(`- Enabled lights: ${scene.lights.filter(l => l.isEnabled()).length}`);
  vrDebug(`- Camera position: ${scene.activeCamera?.position}`);
    
    // Expose debug info globally for Quest console access
    // Provide a teardown to stop VR and dispose HUDs when switching modes
    async function teardownVR() {
      try { engine.stopRenderLoop(); } catch {}
      try { xrHud?.dispose && xrHud.dispose(); } catch {}
      try { vrUI?.dispose && vrUI.dispose(); } catch {}
      try {
        if (vrHelper?.baseExperience?.exitXRAsync) {
          await vrHelper.baseExperience.exitXRAsync();
        } else if (vrHelper?.baseExperience?.sessionManager?.endSessionAsync) {
          await vrHelper.baseExperience.sessionManager.endSessionAsync();
        }
      } catch {}
    }
    // Expose for page-level fallback to desktop
    try { window.vrTeardown = teardownVR; } catch {}

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
        vrDebug("=== VR DEBUG INFO ===");
        vrDebug("Scene meshes:", scene.meshes.map(m => `${m.name} (visible: ${m.isVisible})`));
        vrDebug("Scene lights:", scene.lights.map(l => `${l.constructor.name} (intensity: ${l.intensity}, enabled: ${l.isEnabled()})`));
        vrDebug("Camera:", scene.activeCamera);
        vrDebug("VR Session active:", vrHelper.baseExperience.sessionManager.currentSession !== null);
        vrDebug("Molecule change counter:", mol.changeCounter);
        vrDebug("=== END DEBUG ===");
      }
    };
    
    return {
      engine,
      scene,
      mol,
      state,
      vrHelper,
      vrUI,
      dispose: teardownVR
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
    // Snap positions to exact state by replaying ops from initial to avoid drift; print debug
    try {
      if (window.appState?.recomputeAndCommit) {
        window.appState.recomputeAndCommit();
      }
      if (window.appState?.debugPrint) window.appState.debugPrint('[VR rotate]');
    } catch {}
    // One-shot recompute like desktop: reset flag after use
    if (ui.recompute) ui.recompute = false;
    // Keep bond selected: do not clear window.vrSelection here.
    // Re-emit label update by touching a flag (safe no-op)
    try { scene._lastHudTick = (scene._lastHudTick || 0) + 1; } catch {}
  } catch (e) {
    console.warn('[VR] rotateSelectedBond failed:', e?.message);
  }
}

// createVRInstructions removed (unused)
