// vr/main-vr.js - VR-enabled version of the molecular viewer
import { setupVR } from './vr-setup.js';
import { createVRUI, createVRUIOnCamera } from './vr-ui.js';
// Plot disabled in VR Lite mode

// Import existing modules
import { setupScene, setupMolecule } from '../setup/app-setup.js';
import { createPhysicsManager } from '../physics/physics-manager.js';
import { showMoleculeSelector } from '../ui/molecule-selector.js';

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
    // Step index helper defined immediately after molecule creation so ALL
    // subsequent closures (HUD handlers, XR session callbacks, render loop)
    // capture it. Previous placement caused ReferenceError in some loops.
    function currentStepIndex() {
      try {
        if (typeof mol?.changeCounter === 'number') return mol.changeCounter;
        window.__vrStepSeq = (window.__vrStepSeq || 0) + 1;
        return window.__vrStepSeq;
      } catch {
        return 0;
      }
    }
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
  // Wire molecules selector button (desktop overlay inside VR page)
  try {
    if (vrUI?.btnMolecules) {
      vrUI.btnMolecules.onPointerUpObservable.add(() => {
        try { showMoleculeSelector(); } catch (e) { console.warn('[VR] Failed to open molecule selector', e); }
      });
    }
  } catch {}

  // In-headset molecule panel (XR only) ---------------------------------
  function createXRMoleculePanel(scene, parent) {
    const panelWidth = 0.9; // meters
    const panelHeight = 0.55;
    const root = new BABYLON.TransformNode('xrMolPanelRoot', scene);
    root.parent = parent;
    root.position = new BABYLON.Vector3(0, 0, 0.05); // just slightly in front of HUD root
    const plane = BABYLON.MeshBuilder.CreatePlane('xrMolPanelPlane', { width: panelWidth, height: panelHeight }, scene);
    plane.parent = root;
    plane.isPickable = true;
    const adt = BABYLON.GUI.AdvancedDynamicTexture.CreateForMesh(plane, 2048, 1024, false);
    try { adt.useInvalidateRectOptimization = false; } catch {}
    const bg = new BABYLON.GUI.Rectangle('xrMolPanelBg');
    bg.thickness = 0;
    bg.cornerRadius = 28;
    bg.background = 'rgba(15,18,24,0.92)';
    adt.addControl(bg);
    const title = new BABYLON.GUI.TextBlock('xrMolPanelTitle', 'Select Molecule');
    title.color = '#f5ffd1';
    title.fontSize = 72;
    title.fontWeight = 'bold';
    title.top = '-40%';
    bg.addControl(title);
    const stack = new BABYLON.GUI.StackPanel();
    stack.width = '90%';
    stack.height = '60%';
    stack.top = '5%';
    stack.isVertical = true;
    bg.addControl(stack);
    const status = new BABYLON.GUI.TextBlock('xrMolPanelStatus', 'Loading...');
    status.color = '#cfe3ff';
    status.fontSize = 48;
    status.textWrapping = true;
    stack.addControl(status);
    const btnClose = BABYLON.GUI.Button.CreateSimpleButton('xrMolClose', 'Close');
    btnClose.width = '40%';
    btnClose.height = '14%';
    btnClose.color = '#e9f6ff';
    btnClose.background = 'rgba(60,70,90,0.9)';
    btnClose.cornerRadius = 18;
    btnClose.fontSize = 54;
    btnClose.top = '35%';
    bg.addControl(btnClose);
    btnClose.onPointerUpObservable.add(() => dispose());
    function styleEntry(btn) {
      btn.height = '120px';
      btn.color = '#e9f6ff';
      btn.fontSize = 58;
      btn.thickness = 0;
      btn.background = 'rgba(40,50,65,0.9)';
      btn.cornerRadius = 14;
      btn.onPointerEnterObservable.add(()=> btn.background='rgba(55,70,90,0.95)');
      btn.onPointerOutObservable.add(()=> btn.background='rgba(40,50,65,0.9)');
    }
    async function loadList() {
      try {
        const res = await fetch('/api/molecules');
        if (!res.ok) throw new Error(res.status+'');
        const data = await res.json();
        stack.removeControl(status);
        if (!data.molecules?.length) {
          status.text = 'No .xyz files found';
          stack.addControl(status);
          return;
        }
        data.molecules.forEach(m => {
          const b = BABYLON.GUI.Button.CreateSimpleButton('mol_'+m.name, m.name);
          styleEntry(b);
          b.onPointerUpObservable.add(() => {
            try {
              const base = window.location.pathname;
              const params = new URLSearchParams(window.location.search);
              params.set('molecule', m.name);
              window.location.href = `${base}?${params.toString()}`;
            } catch(e) { console.warn('[VR] reload failed', e); }
          });
          stack.addControl(b);
        });
      } catch (e) {
        status.text = 'Error: '+e.message;
      }
    }
    loadList();
    function dispose() {
      try { adt.dispose(); } catch {}
      try { plane.dispose(); } catch {}
      try { root.dispose(); } catch {}
      window.__xrMolPanel = null;
    }
    return { dispose };
  }

  // Bond rotation controls via 2D HUD overlay
  try {
  const uiState = vrUI?.bond?.state || { side: 'j', step: 5, recompute: false };
  window.vrBondUI = uiState;
  // No bond label in Lite HUD
  let btnMinus = vrUI?.bond?.btnMinus;
  let btnPlus  = vrUI?.bond?.btnPlus;

  // Wire clicks to rotation (with debug)
  btnMinus?.onPointerUpObservable.add(() => { vrDebug('[HUD] (-) pressed'); rotateSelectedBond(scene, -1); try { const e = mlip.compute()?.energy; if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e, currentStepIndex()); } catch {} });
  btnPlus?.onPointerUpObservable.add(() => { vrDebug('[HUD] (+) pressed'); rotateSelectedBond(scene, +1); try { const e = mlip.compute()?.energy; if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e, currentStepIndex()); } catch {} });

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
        // Molecules button (XR)
        try {
          const mBtn = xrHud?.btnMolecules;
          if (mBtn) {
            mBtn.onPointerUpObservable.add(() => {
              if (window.__xrMolPanel) { try { window.__xrMolPanel.dispose(); } catch {}; return; }
              window.__xrMolPanel = createXRMoleculePanel(scene, xrHud.rootMesh || xrCam);
            });
          }
        } catch (e) { console.warn('[VR] Failed wiring XR molecules button', e); }
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
  if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e, currentStepIndex());
      if (xrHud?.energyValue && isFinite(e)) xrHud.energyValue.text = e.toFixed(3);
    } catch {}
  });
  btnPlus?.onPointerUpObservable.add(() => {
    vrDebug('[HUD XR] (+) pressed');
    rotateSelectedBond(scene, +1);
    try {
      const e = mlip.compute()?.energy;
  if (xrHud?.plot && isFinite(e)) xrHud.plot.addPoint(e, currentStepIndex());
      if (xrHud?.energyValue && isFinite(e)) xrHud.energyValue.text = e.toFixed(3);
    } catch {}
  });
  // Seed plot with current energy (if present)
        try {
          const e0 = mlip.compute()?.energy;
          if (xrHud?.plot && isFinite(e0)) xrHud.plot.addPoint(e0, currentStepIndex());
        } catch {}
      } catch {}
      vrLog('[VR] Created camera-anchored HUD for XR session');
    });
    vrHelper.baseExperience.sessionManager.onXRSessionEnded.add(() => {
      try { xrHud?.dispose && xrHud.dispose(); } catch {}
      xrHud = null;
      try { if (window.__xrMolPanel) { window.__xrMolPanel.dispose(); } } catch {}
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
      const physics = createPhysicsManager({ state, getMlip: () => mlip });
      let frameCount = 0;
      let lastEnergyUpdate = 0;
      let lastPlotUpdate = 0;
      const ENERGY_UPDATE_INTERVAL = 5;      // HUD text refresh cadence
      const PLOT_UPDATE_INTERVAL = 20;       // plot point every 20 frames when stable
      let lastPlottedEnergy = null;
      engine.runRenderLoop(() => {
        try {
          frameCount++;
          physics.tickFrame();
          const updated = physics.updatePhysicsCache();
          const e = physics.energy;
          const needsEnergyUpdate = updated || (frameCount - lastEnergyUpdate >= ENERGY_UPDATE_INTERVAL);
          if (needsEnergyUpdate && Number.isFinite(e)) {
            const txt = e.toFixed(3);
            if (vrUI?.energyValue) vrUI.energyValue.text = txt;
            if (xrHud?.energyValue) xrHud.energyValue.text = txt;
            lastEnergyUpdate = frameCount;
          }
          // Plot update: add point on change or periodic interval
          if (xrHud?.plot && Number.isFinite(e)) {
            const energyChanged = (lastPlottedEnergy === null) || Math.abs(e - lastPlottedEnergy) > 1e-9;
            if (updated && energyChanged) {
              xrHud.plot.addPoint(e, currentStepIndex());
              lastPlottedEnergy = e;
              lastPlotUpdate = frameCount;
            } else if ((frameCount - lastPlotUpdate) >= PLOT_UPDATE_INTERVAL) {
              xrHud.plot.addPoint(e, currentStepIndex());
              lastPlottedEnergy = e;
              lastPlotUpdate = frameCount;
            }
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
                  const e2 = physics.energy;
                  if (xrHud?.plot && e2 !== null && isFinite(e2)) {
                    xrHud.plot.addPoint(e2, currentStepIndex());
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
    const r0 = mlip.compute();
    if (r0 && typeof r0.then === 'function') {
      r0.then(res => { if (res && Number.isFinite(res.energy) && vrUI?.energyValue) vrUI.energyValue.text = res.energy.toFixed(3); }).catch(()=>{});
    } else if (r0 && Number.isFinite(r0.energy) && vrUI?.energyValue) {
      vrUI.energyValue.text = r0.energy.toFixed(3);
    }
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
