// vr/vr-setup.js - VR initialization and controller setup
export async function setupVR(engine, scene) {
  console.log('[VR Setup] 🔧 Initializing WebXR...');
  console.log('[VR Setup] Engine ready:', !!engine);
  console.log('[VR Setup] Scene ready:', !!scene);
  console.log('[VR Setup] Scene meshes:', scene.meshes.length);
  
  try {
    // Optimize engine settings for VR performance
    engine.setHardwareScalingLevel(1.0); // Ensure no unnecessary upscaling
    scene.skipPointerMovePicking = true; // Reduce picking overhead
    scene.autoClear = true;
    scene.autoClearDepthAndStencil = true;
    
    console.log('[VR Setup] Engine optimizations applied');
    
    // Try the most minimal approach first - direct WebXR API
    console.log('[VR Setup] Trying direct WebXR API first for better native button support...');
    try {
      const directXR = await BABYLON.WebXRDefaultExperience.CreateAsync(scene, {
        floorMeshes: [],
        optionalFeatures: [], // No optional features
        uiOptions: {
          sessionMode: 'immersive-vr',
          referenceSpaceType: 'local-floor'
        }
      });
      
      if (directXR && directXR.baseExperience) {
        console.log('[VR Setup] ✅ Direct WebXR API successful - this should show native VR button');
        return setupVRFeatures(directXR, scene);
      }
    } catch (directError) {
      console.warn('[VR Setup] Direct WebXR API failed, trying scene method:', directError.message);
    }
    
    // Fallback to scene method
    console.log('[VR Setup] Creating default XR experience via scene...');
  const xrHelper = await scene.createDefaultXRExperienceAsync({
      floorMeshes: [], // No floor mesh needed for molecular viewer
      // Completely disable optional features to prevent errors
      optionalFeatures: [], // No optional features at all
      uiOptions: {
        sessionMode: 'immersive-vr',
        referenceSpaceType: 'local-floor' // Better for Quest headsets
      },
      inputOptions: {
        doNotLoadControllerMeshes: false, // Load controller models
        disableControllerAnimation: false
      },
      // Performance optimizations
      disableTeleportation: true, // We don't need teleportation for molecular viewing
      useStablePlugins: false     // Avoid enabling optional stable plugins automatically
    });
    
    console.log('[VR Setup] XR helper created:', !!xrHelper);
    console.log('[VR Setup] Base experience:', !!xrHelper?.baseExperience);
    console.log('[VR Setup] Session manager:', !!xrHelper?.baseExperience?.sessionManager);
    
    if (!xrHelper.baseExperience) {
      console.warn("[VR Setup] ⚠️ Scene XR initialization failed - trying minimal fallback");
      
      // Try even more minimal initialization
      try {
        console.log('[VR Setup] Attempting absolutely minimal XR setup...');
        const minimalXR = await scene.createDefaultXRExperienceAsync({
          floorMeshes: [],
          optionalFeatures: [],
          disableTeleportation: true,
          uiOptions: {
            sessionMode: 'immersive-vr'
          }
        });
        
        if (minimalXR && minimalXR.baseExperience) {
          console.log('[VR Setup] ✅ Minimal WebXR initialization successful');
          return setupVRFeatures(minimalXR, scene);
        }
      } catch (minimalError) {
        console.error('[VR Setup] Minimal VR initialization failed:', minimalError);
      }
      
      return null;
    }
    
    return setupVRFeatures(xrHelper, scene);
    
  } catch (error) {
    console.error('VR setup failed:', error);
    return null;
  }
}

// Module-scope rotation state used by the polling-based trigger+rotate behavior
let isDragging = false;  // true while a controller trigger is held
let accYaw = 0;          // accumulated yaw rotation (radians)
let accPitch = 0;        // accumulated pitch rotation (radians)

// Helper: get molecule master meshes (atoms and bonds) to rotate directly
function getMoleculeMasters(scene) {
  if (scene._vrMasters && scene._vrMasters.length) return scene._vrMasters;
  const masters = scene.meshes.filter(m => m && m.name && (m.name.startsWith('base_') || m.name.startsWith('bond_')));
  scene._vrMasters = masters;
  console.log('[VR] Found masters to rotate:', masters.map(m => m.name));
  return masters;
}

// Helper: estimate molecule diagonal length for input scaling
function getMoleculeDiag(scene) {
  if (scene._molDiag) return scene._molDiag;
  const masters = getMoleculeMasters(scene);
  if (!masters || masters.length === 0) { scene._molDiag = 1; return 1; }
  let min = new BABYLON.Vector3(Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY);
  let max = new BABYLON.Vector3(Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY);
  for (const m of masters) {
    try {
      m.refreshBoundingInfo && m.refreshBoundingInfo();
      const bi = m.getBoundingInfo();
      const bmin = bi.boundingBox.minimumWorld;
      const bmax = bi.boundingBox.maximumWorld;
      min = BABYLON.Vector3.Minimize(min, bmin);
      max = BABYLON.Vector3.Maximize(max, bmax);
    } catch {}
  }
  const size = max.subtract(min);
  const diag = Math.max(size.x, size.y, size.z) || 1;
  scene._molDiag = diag;
  return diag;
}

function setupVRFeatures(xrHelper, scene) {
  console.log('[VR Setup] Configuring VR features...');
  // Cache masters list for rotation-based controls
  try { getMoleculeMasters(scene); } catch (e) { console.warn('[VR] Could not collect molecule masters:', e?.message); }
  // Track last two-hand center for dual-trigger translation
  let lastTwoHandCenter = null;
  
  // Disable automatic locomotion features - we'll use manual trigger + rotate behavior
  if (xrHelper.baseExperience.featuresManager) {
    console.log('[VR Setup] Skipping optional locomotion/plugins; using simple trigger+rotate molecule control');
    
  } else {
    console.warn('[VR Setup] ⚠️ Features manager not available - using basic VR only');
  }
  
  // VR-specific settings
  xrHelper.baseExperience.sessionManager.onXRSessionInit.add(() => {
    // Reset rotation/behavior flags on session start
    isDragging = false;
    accYaw = 0;
    accPitch = 0;
    scene._behaviorsActive = false;
    scene._grabActive = false;
    console.log("VR Session started successfully");
    
    // Enhanced debugging for black screen issues
    console.log('[VR Debug] Scene state at VR start:');
    console.log('- Meshes:', scene.meshes.length);
    console.log('- Lights:', scene.lights.length);
    console.log('- Active Camera:', scene.activeCamera?.constructor.name);
    console.log('- Camera Position:', scene.activeCamera?.position);
    
    // Check and fix common lighting issues
    let lightIntensityTotal = 0;
    scene.lights.forEach((light, i) => {
      console.log(`- Light ${i}: ${light.constructor.name}, intensity: ${light.intensity}, enabled: ${light.isEnabled()}`);
      lightIntensityTotal += light.intensity;
      
      // Ensure lights are enabled and have adequate intensity for VR
      if (light instanceof BABYLON.HemisphericLight) {
        light.intensity = Math.max(light.intensity, 0.8); // Minimum intensity for VR
        light.setEnabled(true);
      } else if (light instanceof BABYLON.DirectionalLight) {
        light.intensity = Math.max(light.intensity, 1.0);
        light.setEnabled(true);
      }
    });
    
    console.log(`Total light intensity: ${lightIntensityTotal}`);
    
    if (lightIntensityTotal < 0.5) {
      console.warn('[VR Debug] WARNING: Very low total light intensity, adding emergency lighting');
      
      // Add emergency lighting if scene is too dark
      const emergencyLight = new BABYLON.HemisphericLight("emergencyVRLight", new BABYLON.Vector3(0, 1, 0), scene);
      emergencyLight.intensity = 2.0; // Brighter emergency lighting
      emergencyLight.diffuse = new BABYLON.Color3(1, 1, 1);
      console.log('[VR Debug] Added emergency hemisphere light with intensity 2.0');
      
      // Also add a directional emergency light
      const emergencyDir = new BABYLON.DirectionalLight("emergencyVRDir", new BABYLON.Vector3(1, -1, 1), scene);
      emergencyDir.intensity = 1.5;
      emergencyDir.diffuse = new BABYLON.Color3(1, 1, 1);
      console.log('[VR Debug] Added emergency directional light');
    } else {
      // Even if we have lights, boost them for VR
      console.log('[VR Debug] Boosting existing lights for VR visibility');
      scene.lights.forEach((light, i) => {
        const originalIntensity = light.intensity;
        if (light instanceof BABYLON.HemisphericLight) {
          light.intensity = Math.max(light.intensity, 1.5);
        } else if (light instanceof BABYLON.DirectionalLight) {
          light.intensity = Math.max(light.intensity, 1.2);
        }
        if (light.intensity !== originalIntensity) {
          console.log(`[VR Debug] Boosted light ${i} from ${originalIntensity} to ${light.intensity}`);
        }
      });
    }
    
    // Check camera positioning for VR
    if (scene.activeCamera) {
      const camera = scene.activeCamera;
      console.log('[VR Debug] Camera details:');
      console.log('- Type:', camera.constructor.name);
      console.log('- Position:', camera.position);
      console.log('- Target:', camera.target || 'N/A');
      console.log('- FOV:', camera.fov || 'N/A');
      
      // Ensure camera is positioned appropriately for molecules
      if (camera.position.length() < 0.1) {
        console.warn('[VR Debug] Camera too close to origin, adjusting...');
        camera.position = new BABYLON.Vector3(0, 0, 5);
      }
      
      // Set target to molecular center
      if (camera.setTarget) {
        camera.setTarget(BABYLON.Vector3.Zero());
      }
    }
    
    // Check mesh visibility and materials
    let visibleMeshCount = 0;
    scene.meshes.forEach((mesh, i) => {
      if (mesh.isVisible && mesh.isEnabled()) {
        visibleMeshCount++;
        
        // Check material
        if (!mesh.material) {
          console.warn(`[VR Debug] Mesh ${i} (${mesh.name}) has no material - adding default`);
          // Add default material
          const defaultMaterial = new BABYLON.StandardMaterial(`defaultVR_${mesh.name}`, scene);
          defaultMaterial.diffuseColor = new BABYLON.Color3(0.8, 0.8, 0.8);
          defaultMaterial.emissiveColor = new BABYLON.Color3(0.1, 0.1, 0.1); // Self-illuminated
          mesh.material = defaultMaterial;
        } else {
          // Ensure material is visible in VR
          if (mesh.material.alpha < 0.1) {
            console.warn(`[VR Debug] Mesh ${i} (${mesh.name}) is nearly transparent (alpha: ${mesh.material.alpha}) - fixing`);
            mesh.material.alpha = 1.0;
          }
          
          // Add some emissive lighting to materials to ensure they're visible even with poor lighting
          if (mesh.material instanceof BABYLON.StandardMaterial && mesh.material.emissiveColor) {
            const originalEmissive = mesh.material.emissiveColor.clone();
            mesh.material.emissiveColor = new BABYLON.Color3(
              Math.max(originalEmissive.r, 0.05),
              Math.max(originalEmissive.g, 0.05),
              Math.max(originalEmissive.b, 0.05)
            );
          }
        }
      }
    });
    
    console.log(`[VR Debug] Visible meshes: ${visibleMeshCount}/${scene.meshes.length}`);
    
    if (visibleMeshCount === 0) {
      console.error('[VR Debug] NO VISIBLE MESHES - This will cause a black screen!');
      
      // Create a test cube to verify rendering is working
      const testCube = BABYLON.MeshBuilder.CreateBox("vrTestCube", {size: 1}, scene);
      testCube.position = new BABYLON.Vector3(0, 0, -3);
      const testMaterial = new BABYLON.StandardMaterial("vrTestMaterial", scene);
      testMaterial.diffuseColor = new BABYLON.Color3(1, 0, 0); // Red cube
      testMaterial.emissiveColor = new BABYLON.Color3(0.2, 0, 0); // Self-illuminated
      testCube.material = testMaterial;
      console.log('[VR Debug] Created red test cube at (0, 0, -3)');
    }
    
    // Force render to ensure everything is updated
    scene.render();
    console.log('[VR Debug] Forced initial render');
    
    // Add continuous lighting monitoring to prevent lights being turned off
    let lightingCheckInterval = setInterval(() => {
      let currentLightTotal = 0;
      let enabledLights = 0;
      scene.lights.forEach(light => {
        if (light.isEnabled()) {
          currentLightTotal += light.intensity;
          enabledLights++;
        }
      });
      
      if (currentLightTotal < 1.0 || enabledLights === 0) {
        console.warn(`[VR Debug] Lighting dimmed/disabled! Total: ${currentLightTotal}, Enabled: ${enabledLights} - Restoring...`);
        
        // Re-enable and boost all lights
        scene.lights.forEach(light => {
          light.setEnabled(true);
          if (light instanceof BABYLON.HemisphericLight) {
            light.intensity = Math.max(light.intensity, 1.5);
          } else if (light instanceof BABYLON.DirectionalLight) {
            light.intensity = Math.max(light.intensity, 1.2);
          }
        });
        
        // Add emergency light if needed
        if (enabledLights === 0) {
          const emergencyRestore = new BABYLON.HemisphericLight("emergencyRestore_" + Date.now(), new BABYLON.Vector3(0, 1, 0), scene);
          emergencyRestore.intensity = 2.0;
          emergencyRestore.diffuse = new BABYLON.Color3(1, 1, 1);
          console.log('[VR Debug] Added emergency restore light');
        }
      }
    }, 2000); // Check every 2 seconds
    
    // Stop monitoring when VR session ends
    xrHelper.baseExperience.sessionManager.onXRSessionEnded.addOnce(() => {
      clearInterval(lightingCheckInterval);
      console.log('[VR Debug] Stopped lighting monitoring');
    });
    
    // Show the browser's native VR button by triggering WebXR state
    console.log('[VR Debug] VR session active - native VR controls should be available');

    // Try enabling engine-provided grab+scale behaviors ("standard" interaction)
    try {
      if (window && window.vrDisableBehaviors) {
        console.log('[VR Setup] XR behaviors disabled via flag');
      }
      const activated = !window?.vrDisableBehaviors && enableStandardXRInteractions(scene, xrHelper);
      if (activated) {
        scene._behaviorsActive = true;
        console.log('[VR Setup] ✅ Standard XR behaviors active: Grab to rotate/move, two hands to scale');
      } else {
        console.log('[VR Setup] ℹ️ Standard XR behaviors not activated; falling back to trigger-rotate');
      }
    } catch (bx) {
      console.warn('[VR Setup] XR behaviors setup failed; using fallback trigger-rotate:', bx?.message);
      scene._behaviorsActive = false;
    }
  });
  
  xrHelper.baseExperience.sessionManager.onXRSessionEnded.add(() => {
    console.log("VR Session ended");
    // Cleanup behavior wrapper and flags
    try {
      if (scene._molWrapper && scene._molWrapper.dispose) {
        scene._molWrapper.dispose();
      }
      if (scene._behaviorObserver) {
        scene.onBeforeRenderObservable.remove(scene._behaviorObserver);
      }
      delete scene._molWrapper;
      delete scene._behaviorObserver;
      delete scene._molRotation;
      delete scene._molScale;
      window.vrBehaviorsActive = false;
      scene._behaviorsActive = false;
      scene._grabActive = false;
    } catch {}
  });
  
  // Per-controller polling fallback (works even when observables/components are missing)
  const controllerState = new WeakMap(); // controller -> { pressed, lastYaw, lastPitch }
  const getControllerNode = (controller) => controller.pointer || controller.grip || controller.pointer?.parent || null;
  const getYawPitchForController = (controller, xrHelper) => {
    const node = getControllerNode(controller);
    if (node && node.getDirection) {
      const f = node.getDirection(BABYLON.Vector3.Forward()).normalize();
      return { yaw: Math.atan2(f.x, f.z), pitch: Math.asin(-Math.max(-1, Math.min(1, f.y))) };
    }
    // Fallback: WebXR pose from targetRaySpace
    try {
      const inputSource = controller.inputSource;
      const sessionMan = xrHelper.baseExperience.sessionManager;
      const frame = sessionMan.currentFrame;
      const ref = sessionMan.referenceSpace;
      if (inputSource?.targetRaySpace && frame && ref) {
        const pose = frame.getPose(inputSource.targetRaySpace, ref);
        if (pose?.transform?.matrix) {
          const m = pose.transform.matrix;
          const dir = new BABYLON.Vector3(-m[8], -m[9], -m[10]); // -Z is forward
          const f = dir.normalize();
          return { yaw: Math.atan2(f.x, f.z), pitch: Math.asin(-Math.max(-1, Math.min(1, f.y))) };
        }
      }
    } catch {}
    return { yaw: 0, pitch: 0 };
  };

  const getControllerWorldPosition = (controller, xrHelper) => {
    try {
      const node = getControllerNode(controller);
      if (node && node.getAbsolutePosition) {
        return node.getAbsolutePosition().clone();
      }
      const inputSource = controller.inputSource;
      const sessionMan = xrHelper.baseExperience.sessionManager;
      const frame = sessionMan.currentFrame;
      const ref = sessionMan.referenceSpace;
      if (inputSource?.targetRaySpace && frame && ref) {
        const pose = frame.getPose(inputSource.targetRaySpace, ref);
        if (pose?.transform) {
          if (pose.transform.position) {
            return new BABYLON.Vector3(
              pose.transform.position.x,
              pose.transform.position.y,
              pose.transform.position.z
            );
          }
          const m = pose.transform.matrix;
          if (m && m.length >= 16) {
            return new BABYLON.Vector3(m[12], m[13], m[14]);
          }
        }
      }
    } catch {}
    return null;
  };

  const isTriggerPressed = (controller) => {
    try {
      const mc = controller.motionController || null;
      if (mc && mc.getComponent) {
        const trig = mc.getComponent('xr-standard-trigger');
        if (trig && typeof trig.pressed === 'boolean') return trig.pressed;
      }
      const gp = controller.inputSource?.gamepad;
      if (gp && gp.buttons && gp.buttons.length > 0) {
        // Heuristic: button[0] is trigger on many controllers
        return !!gp.buttons[0]?.pressed;
      }
    } catch {}
    return false;
  };

  xrHelper.input.onControllerAddedObservable.add((controller) => {
    console.log('VR Controller added:', controller.inputSource.handedness, controller.inputSource.profiles);
    controllerState.set(controller, { pressed: false, lastYaw: 0, lastPitch: 0 });
  });

  // Scene-level polling loop
  scene.onBeforeRenderObservable.add(() => {
    // Skip polling-based rotation only while behaviors are actively dragging/scaling
    if (scene._behaviorsActive && scene._grabActive) return;
    const masters = getMoleculeMasters(scene);
    if (!masters || masters.length === 0) return;
    const inputs = xrHelper.input?.controllers || [];

    // Determine if both triggers are held for two-hand translation
    const pressedWithPos = [];
    for (const c of inputs) {
      if (isTriggerPressed(c)) {
        const p = getControllerWorldPosition(c, xrHelper);
        if (p) pressedWithPos.push({ c, p });
      }
    }
    const twoHandMode = pressedWithPos.length >= 2;
    for (const controller of inputs) {
      let st = controllerState.get(controller);
      if (!st) { st = { pressed: false, lastYaw: 0, lastPitch: 0 }; controllerState.set(controller, st); }
      const pressed = isTriggerPressed(controller);
      if (pressed && !st.pressed) {
        // Edge: press start
        const { yaw, pitch } = getYawPitchForController(controller, xrHelper);
        st.lastYaw = yaw; st.lastPitch = pitch; st.pressed = true;
        // mark global dragging so we accumulate deltas
        isDragging = true;
        // console.log('[VR] trigger down', controller.inputSource.handedness || 'unknown');
      } else if (!pressed && st.pressed) {
        // Edge: press end
        st.pressed = false;
        isDragging = false;
        // console.log('[VR] trigger up');
      }
      if (st.pressed && !twoHandMode) {
        const sensitivity = 0.9;
        const { yaw, pitch } = getYawPitchForController(controller, xrHelper);
        let dYaw = yaw - st.lastYaw;
        let dPitch = pitch - st.lastPitch;
        if (dYaw > Math.PI) dYaw -= 2 * Math.PI; else if (dYaw < -Math.PI) dYaw += 2 * Math.PI;
        if (dPitch > Math.PI) dPitch -= 2 * Math.PI; else if (dPitch < -Math.PI) dPitch += 2 * Math.PI;
        const maxStep = 0.05;
        dYaw = Math.max(-maxStep, Math.min(maxStep, dYaw));
        dPitch = Math.max(-maxStep, Math.min(maxStep, dPitch));
        accYaw += -dYaw * sensitivity;
        accPitch += -dPitch * sensitivity;
        const q = BABYLON.Quaternion.RotationAxis(BABYLON.Axis.Y, accYaw).multiply(
                  BABYLON.Quaternion.RotationAxis(BABYLON.Axis.X, accPitch));
        for (const m of masters) m.rotationQuaternion = q.clone();
        st.lastYaw = yaw; st.lastPitch = pitch;
      }
    }

    // Apply two-hand translation based on average controller position delta
    if (twoHandMode) {
      const center = pressedWithPos[0].p.add(pressedWithPos[1].p).scale(0.5);
      if (!lastTwoHandCenter) {
        lastTwoHandCenter = center.clone();
      } else {
        // Compute delta and amplify sensitivity
        let delta = center.subtract(lastTwoHandCenter);
        const cam = scene.activeCamera;
        // Optionally convert to head-relative space for intuitive push/pull
        if (cam && cam.getDirection) {
          const right = cam.getDirection(BABYLON.Axis.X);
          const up = cam.getDirection(BABYLON.Axis.Y);
          const fwd = cam.getDirection(BABYLON.Axis.Z);
          // Project delta onto camera axes and scale
          const dx = BABYLON.Vector3.Dot(delta, right);
          const dy = BABYLON.Vector3.Dot(delta, up);
          const dz = BABYLON.Vector3.Dot(delta, fwd);
          // Sensitivity gains (allow runtime override via window.vrTwoHandGains)
          const override = (typeof window !== 'undefined' && window.vrTwoHandGains) ? window.vrTwoHandGains : null;
          const gainF = override?.f ?? 60.0; // push/pull (5x prior 12)
          const gainL = override?.l ?? 16.0; // lateral (2x prior 8)
          const gainU = override?.u ?? 16.0; // vertical (2x prior 8)
          delta = right.scale(dx * gainL).addInPlace(up.scale(dy * gainU)).addInPlace(fwd.scale(dz * gainF));
        } else {
          delta.scaleInPlace(5.0); // generic gain if camera dirs not available
        }
        // Relax clamp to allow bigger moves while keeping stability
        const clampOverride = (typeof window !== 'undefined' && window.vrTwoHandClampFactor) ? window.vrTwoHandClampFactor : null;
        const maxStep = getMoleculeDiag(scene) * (clampOverride ?? 2.0);
        if (delta.length() > maxStep) {
          delta.normalize().scaleInPlace(maxStep);
        }
        if (scene._applyMolTranslation) {
          scene._applyMolTranslation(delta);
        } else {
          for (const m of masters) m.position.addInPlace(delta);
        }
        lastTwoHandCenter.copyFrom(center);
      }
    } else {
      lastTwoHandCenter = null;
    }

    // Joystick pan and zoom (when not grabbing and not pressing trigger)
    try {
      const deadzone = 0.1;
      const diag = getMoleculeDiag(scene);
  const panSpeed = Math.max(0.01, 0.03) * diag; // meters per full deflection per frame
  const zoomGain = 0.125; // half the previous sensitivity
      let leftX = 0, leftY = 0, rightY = 0;
      for (const controller of inputs) {
        const pressed = isTriggerPressed(controller);
        if (pressed) continue; // don't mix with trigger actions
        const gp = controller.inputSource?.gamepad;
        if (!gp || !gp.axes || gp.axes.length < 2) continue;
        const ax = gp.axes[0] || 0;
        const ay = gp.axes[1] || 0;
        const hand = controller.inputSource?.handedness;
        if (hand === 'left') {
          leftX = Math.abs(ax) > deadzone ? ax : 0;
          leftY = Math.abs(ay) > deadzone ? ay : 0;
        } else if (hand === 'right') {
          // Prefer right stick Y if present (indices 3), else fallback to 1
          const ayr = gp.axes.length >= 4 ? (gp.axes[3] || 0) : ay;
          rightY = Math.abs(ayr) > deadzone ? ayr : 0;
        }
      }

      // Pan: left stick X (camera right), Y (camera up)
      if (leftX !== 0 || leftY !== 0) {
        const cam = scene.activeCamera;
        const right = cam && cam.getDirection ? cam.getDirection(BABYLON.Axis.X) : BABYLON.Axis.X;
        const up = cam && cam.getDirection ? cam.getDirection(BABYLON.Axis.Y) : BABYLON.Axis.Y;
        // Invert X so pushing left moves left in world
        const delta = right.scale(leftX * panSpeed).addInPlace(up.scale(leftY * panSpeed));
        if (scene._applyMolTranslation) {
          scene._applyMolTranslation(delta);
        } else {
          for (const m of masters) m.position.addInPlace(delta);
        }
      }

      // Zoom: right stick Y (scale molecule). Up to zoom-in (reduce ay negative => factor > 1)
      if (rightY !== 0) {
  const factor = Math.exp(-rightY * zoomGain);
  const clamped = Math.max(0.95, Math.min(1.05, factor));
        if (scene._applyMolScaleFactor) {
          scene._applyMolScaleFactor(clamped);
        } else {
          for (const m of masters) {
            const s = m.scaling || new BABYLON.Vector3(1,1,1);
            m.scaling = new BABYLON.Vector3(s.x * clamped, s.y * clamped, s.z * clamped);
          }
        }
      }
    } catch {}
  });
  
  // Add error handling for missing features - but don't try to enable hand tracking
  if (xrHelper.baseExperience.featuresManager) {
    console.log('[VR Setup] Features manager available');
    
    // Skip hand tracking to avoid initialization errors
    console.log('[VR Setup] Skipping optional features to ensure VR button appears');
  }
  
  // CRITICAL: Properly set up WebXR for native browser button
  setTimeout(() => {
    console.log('[VR Setup] Triggering native VR button availability...');
    
    try {
      // Method 1: Ensure the XR helper is properly attached to the scene
      if (xrHelper.baseExperience && xrHelper.baseExperience.sessionManager) {
        console.log('[VR Setup] Session manager properly initialized');
        
        // Method 2: Check if browser can detect VR capability
        if (navigator.xr) {
          navigator.xr.isSessionSupported('immersive-vr').then(supported => {
            if (supported) {
              console.log('[VR Setup] ✅ Immersive VR supported - native Enter VR button should appear');
              
              // Method 3: Just ensure user activation without creating sessions
              console.log('[VR Setup] Ensuring user activation for native VR button...');
              
              // The WebXR system is ready, browser should show native button
              // No need to create temporary sessions that enter VR mode
              const canvas = scene.getEngine().getRenderingCanvas();
              if (canvas) {
                // Add one-time click listener to ensure user activation
                const activateVR = (event) => {
                  console.log('[VR Setup] ✅ User interaction detected - native VR button should be active');
                  canvas.removeEventListener('click', activateVR);
                  canvas.removeEventListener('touchstart', activateVR);
                };
                
                canvas.addEventListener('click', activateVR, { once: true });
                canvas.addEventListener('touchstart', activateVR, { once: true });
                
                console.log('[VR Setup] 🎯 User activation listeners added - VR button should appear after any interaction');
              }
              
              // Also dispatch a custom event to signal VR readiness
              window.dispatchEvent(new CustomEvent('webxr-immersive-ready', {
                detail: { xrHelper, canvas }
              }));
              
            } else {
              console.warn('[VR Setup] ⚠️ Immersive VR not supported');
            }
          }).catch(err => {
            console.error('[VR Setup] VR support check failed:', err);
          });
        }
      }
    } catch (error) {
      console.error('[VR Setup] Native button trigger failed:', error);
    }
  }, 2000); // Give WebXR more time to fully initialize
  
  return xrHelper;
}

// Export the main setup function
export { setupVRFeatures };

// ===========
// XR Behaviors (optional): Grab to rotate, two-hand pinch to scale
// ===========
function enableStandardXRInteractions(scene, xrHelper) {
  try {
    const fm = xrHelper?.baseExperience?.featuresManager;
    if (!fm) {
      console.warn('[VR Behaviors] No featuresManager; cannot enable pointer selection');
      return false;
    }
    // Try to enable pointer selection for controller rays
    try {
      fm.enableFeature(BABYLON.WebXRFeatureName.POINTER_SELECTION, 'latest', {
        xrInput: xrHelper.input,
        forceGazeMode: false,
        disableScenePointerVectorUpdate: false
      });
      // Near interaction can improve grabbing with close controllers (optional)
      try { fm.enableFeature(BABYLON.WebXRFeatureName.NEAR_INTERACTION, 'latest'); } catch {}
    } catch (pfErr) {
      console.warn('[VR Behaviors] Pointer selection enable failed:', pfErr?.message);
      return false;
    }

    // Build an invisible wrapper box around the molecule to act as a grab target
    const masters = getMoleculeMasters(scene);
    if (!masters || masters.length === 0) return false;

    // Compute molecule bounds in world space
    let min = new BABYLON.Vector3(Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY, Number.POSITIVE_INFINITY);
    let max = new BABYLON.Vector3(Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY, Number.NEGATIVE_INFINITY);
    for (const m of masters) {
      try {
        m.refreshBoundingInfo && m.refreshBoundingInfo();
        const bi = m.getBoundingInfo();
        const bmin = bi.boundingBox.minimumWorld;
        const bmax = bi.boundingBox.maximumWorld;
        min = BABYLON.Vector3.Minimize(min, bmin);
        max = BABYLON.Vector3.Maximize(max, bmax);
      } catch {}
    }
    const center = min.add(max).scale(0.5);
    const size = max.subtract(min);
    const diag = Math.max(size.x, size.y, size.z) || 1;

    const inflate = 1.1; // enlarge slightly to make grabbing easier
    const wrapper = BABYLON.MeshBuilder.CreateBox('molWrapper', {
      width: (size.x || 0.1) * inflate,
      height: (size.y || 0.1) * inflate,
      depth: (size.z || 0.1) * inflate
    }, scene);
    wrapper.position.copyFrom(center);
    wrapper.rotationQuaternion = BABYLON.Quaternion.Identity();
    wrapper.isPickable = true;
    try { wrapper.isNearPickable = true; } catch {}
    try { wrapper.enablePointerMoveEvents = true; } catch {}
    wrapper.visibility = 0.02; // almost invisible but grabbable; set to 0 to hide completely
    wrapper.alwaysSelectAsActiveMesh = true;
    wrapper.checkCollisions = false;

    // Add grab (6DoF) and two-hand scale behaviors
    const sixDof = new BABYLON.SixDofDragBehavior();
    try {
      sixDof.onDragStartObservable.add(() => { scene._grabActive = true; console.log('[VR Behaviors] Drag start'); });
      sixDof.onDragObservable.add(() => { scene._grabActive = true; });
      sixDof.onDragEndObservable.add(() => { scene._grabActive = false; console.log('[VR Behaviors] Drag end'); });
    } catch {}
    wrapper.addBehavior(sixDof);
    const scaler = new BABYLON.MultiPointerScaleBehavior();
    wrapper.addBehavior(scaler);

    // Accumulated molecule transform
  let molRotation = BABYLON.Quaternion.Identity();
  let molScale = 1;
  let currentCenter = center.clone();
  let prevWrapperPos = wrapper.position.clone();

    // Initialize masters with rotationQuaternion
    for (const m of masters) {
      if (!m.rotationQuaternion) m.rotationQuaternion = BABYLON.Quaternion.Identity();
      m.scaling = new BABYLON.Vector3(molScale, molScale, molScale);
    }

    // Use a per-frame poll of wrapper deltas to update masters, then reset wrapper
    const epsilon = 1e-6;
    const behaviorObserver = scene.onBeforeRenderObservable.add(() => {
      let activeThisFrame = false;
      // Rotation delta
      if (wrapper.rotationQuaternion) {
        const rq = wrapper.rotationQuaternion;
        if (Math.abs(1 - rq.length()) > epsilon || Math.abs(rq.x) > epsilon || Math.abs(rq.y) > epsilon || Math.abs(rq.z) > epsilon) {
          // Apply delta then reset wrapper rotation
          molRotation = rq.multiply(molRotation);
          for (const m of masters) m.rotationQuaternion = molRotation.clone();
          wrapper.rotationQuaternion = BABYLON.Quaternion.Identity();
          activeThisFrame = true;
        }
      }
      // Translation delta
      if (wrapper.position && prevWrapperPos) {
        const delta = wrapper.position.subtract(prevWrapperPos);
        if (delta.lengthSquared() > 1e-6) {
          currentCenter.addInPlace(delta);
          for (const m of masters) m.position.addInPlace(delta);
          prevWrapperPos.copyFrom(wrapper.position);
          activeThisFrame = true;
        }
        // Keep wrapper centered on molecule after applying
        if (!wrapper.position.equalsWithEpsilon(currentCenter, 1e-6)) {
          wrapper.position.copyFrom(currentCenter);
          prevWrapperPos.copyFrom(currentCenter);
        }
      }
      // Scale delta (use X as representative; assume uniform)
      const s = wrapper.scaling;
      if (s && (Math.abs(s.x - 1) > 0.001 || Math.abs(s.y - 1) > 0.001 || Math.abs(s.z - 1) > 0.001)) {
        const scaleFactor = Math.max(0.25, Math.min(4, molScale * s.x));
        molScale = scaleFactor;
        const sv = new BABYLON.Vector3(molScale, molScale, molScale);
        for (const m of masters) m.scaling = sv;
        wrapper.scaling.copyFromFloats(1, 1, 1);
        activeThisFrame = true;
      }
      // Expose whether behaviors were actively manipulating this frame
      scene._grabActive = !!activeThisFrame;
  });

    // Helpers to apply transforms consistently from any input source
    scene._applyMolTranslation = (delta) => {
      currentCenter.addInPlace(delta);
      for (const m of masters) m.position.addInPlace(delta);
      wrapper.position.copyFrom(currentCenter);
      prevWrapperPos.copyFrom(currentCenter);
    };
    scene._applyMolScaleFactor = (factor) => {
      molScale = Math.max(0.25, Math.min(4, molScale * factor));
      const sv = new BABYLON.Vector3(molScale, molScale, molScale);
      for (const m of masters) m.scaling = sv;
    };

    // Store on scene for debug
  scene._molWrapper = wrapper;
    scene._molRotation = () => molRotation.clone();
    scene._molScale = () => molScale;
  scene._behaviorObserver = behaviorObserver;
    window.vrBehaviorsActive = true;
    return true;
  } catch (err) {
    console.warn('[VR Behaviors] Setup failed:', err?.message || err);
    window.vrBehaviorsActive = false;
    return false;
  }
}
