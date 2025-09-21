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
      useStablePlugins: true      // Use stable, tested plugins only
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

function setupVRFeatures(xrHelper, scene) {
  console.log('[VR Setup] Configuring VR features...');
  
  // Disable automatic locomotion features - we'll use manual trigger + drag camera control
  if (xrHelper.baseExperience.featuresManager) {
    console.log('[VR Setup] Skipping automatic locomotion - using manual trigger + drag camera control...');
    
    // Note: We're intentionally NOT enabling teleportation or movement features
    // to avoid conflicts with our manual camera control system
    
    try {
      // Add pointer selection for better bond interaction
      const pointerFeature = xrHelper.baseExperience.featuresManager.enableFeature(BABYLON.WebXRFeatureName.POINTER_SELECTION, "stable", {
        xrInput: xrHelper.input,
        enablePointerSelectionOnAllControllers: true,
        preferredHandedness: "none", // Enable on both hands
        disablePointerUpOnTouchDown: false
      });
      
      if (pointerFeature) {
        console.log('[VR Setup] ✅ Pointer selection enabled - point and click to select bonds');
        
        // Set up pointer selection events
        pointerFeature.onPointerPickObservable.add((pickingInfo) => {
          if (pickingInfo.hit && pickingInfo.pickedMesh) {
            console.log('[VR] 🎯 Pointer picked:', pickingInfo.pickedMesh.name);
            
            // Check if it's a bond
            if (pickingInfo.pickedMesh.name && pickingInfo.pickedMesh.name.startsWith("bond_")) {
              handleBondSelection(pickingInfo.pickedMesh);
            }
          }
        });
        
        // Also add selection events
        pointerFeature.onPointerSelectionObservable.add((selectionInfo) => {
          if (selectionInfo.hit && selectionInfo.pickedMesh) {
            console.log('[VR] 🎯 Selection event:', selectionInfo.pickedMesh.name);
            
            // Check if it's a bond
            if (selectionInfo.pickedMesh.name && selectionInfo.pickedMesh.name.startsWith("bond_")) {
              handleBondSelection(selectionInfo.pickedMesh);
            }
          }
        });
      }
    } catch (pointerError) {
      console.warn('[VR Setup] ⚠️ Pointer selection feature failed:', pointerError.message);
      console.log('[VR Setup] Will use manual controller interaction as fallback');
    }
    
    // Add physics controller interaction as backup
    try {
      const physicsFeature = xrHelper.baseExperience.featuresManager.enableFeature(BABYLON.WebXRFeatureName.PHYSICS_CONTROLLERS, "stable", {
        xrInput: xrHelper.input,
        physicsProperties: {
          restitution: 0.5,
          impostorSize: 0.1
        }
      });
      
      if (physicsFeature) {
        console.log('[VR Setup] ✅ Physics controllers enabled for interaction');
      }
    } catch (physicsError) {
      console.warn('[VR Setup] ⚠️ Physics controllers failed:', physicsError.message);
    }
    
  } else {
    console.warn('[VR Setup] ⚠️ Features manager not available - using basic VR only');
  }
  
  // VR-specific settings
  xrHelper.baseExperience.sessionManager.onXRSessionInit.add(() => {
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
  });
  
  xrHelper.baseExperience.sessionManager.onXRSessionEnded.add(() => {
    console.log("VR Session ended");
  });
  
  // Enhanced hand controller setup for Meta Quest with error handling
  xrHelper.input.onControllerAddedObservable.add((controller) => {
    console.log("VR Controller added:", controller.inputSource.handedness, controller.inputSource.profiles);
    
    controller.onMotionControllerInitObservable.add((motionController) => {
      console.log('Motion controller initialized:', motionController.profile);
      try {
        setupControllerInteraction(motionController, scene);
      } catch (controllerError) {
        console.warn('Controller setup failed, but VR will still work:', controllerError);
      }
    });
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

// Global camera drag state (shared across controller functions)
let isDragging = false;
let lastControllerPosition = null;
let lastControllerRotation = null;

function setupControllerInteraction(motionController, scene) {
  console.log('[VR] Setting up controller interaction for:', motionController.profile);
  console.log('[VR] Controller properties:', {
    rootMesh: !!motionController.rootMesh,
    mesh: !!motionController.mesh,
    inputSource: !!motionController.inputSource,
    hasComponents: !!motionController.getComponent
  });
  
  try {
    // Get the main trigger component
    const triggerComponent = motionController.getComponent("xr-standard-trigger");
    const squeezeComponent = motionController.getComponent("xr-standard-squeeze");
    
    if (triggerComponent) {
      console.log('[VR] ✅ Trigger component found - will use for camera control + bond selection');
      
      // Trigger for camera control (drag) and bond selection (quick press)
      let triggerPressTime = 0;
      let triggerStartPosition = null;
      
      triggerComponent.onButtonStateChangedObservable.add((component) => {
        if (component.pressed) {
          triggerPressTime = Date.now();
          isDragging = true;
          
          // Get initial controller position for drag calculation
          const controllerMesh = motionController.rootMesh || motionController.mesh;
          if (controllerMesh) {
            triggerStartPosition = controllerMesh.absolutePosition?.clone() || controllerMesh.position?.clone();
            lastControllerPosition = triggerStartPosition?.clone();
            lastControllerRotation = controllerMesh.rotationQuaternion?.clone() || 
              BABYLON.Quaternion.FromEulerAngles(controllerMesh.rotation.x, controllerMesh.rotation.y, controllerMesh.rotation.z);
          }
          
          console.log('[VR] 🎯 Trigger pressed - starting camera drag or bond selection');
          
        } else {
          // Trigger released
          const pressDuration = Date.now() - triggerPressTime;
          const controllerMesh = motionController.rootMesh || motionController.mesh;
          
          isDragging = false;
          
          // If it was a quick press (< 300ms) and minimal movement, treat as bond selection
          if (pressDuration < 300) {
            let totalMovement = 0;
            if (triggerStartPosition && controllerMesh) {
              const currentPos = controllerMesh.absolutePosition || controllerMesh.position;
              totalMovement = BABYLON.Vector3.Distance(triggerStartPosition, currentPos);
            }
            
            if (totalMovement < 0.1) { // Less than 10cm movement = selection
              console.log('[VR] 🎯 Quick trigger press - attempting bond selection');
              performBondSelection(motionController, scene);
            } else {
              console.log('[VR] 📹 Trigger drag completed - camera moved');
            }
          } else {
            console.log('[VR] 📹 Long trigger press - camera drag completed');
          }
          
          // Reset tracking variables
          triggerStartPosition = null;
          lastControllerPosition = null;
          lastControllerRotation = null;
        }
      });
    } else {
      console.warn('[VR] ⚠️ No trigger component found on controller');
    }
    
    if (squeezeComponent) {
      console.log('[VR] ✅ Squeeze component found - will use for menu');
      // Squeeze for menu activation
      squeezeComponent.onButtonStateChangedObservable.add((component) => {
        if (component.pressed) {
          console.log('[VR] 🎯 Squeeze pressed - toggling menu');
          toggleVRMenu(motionController, scene);
        }
      });
    } else {
      console.warn('[VR] ⚠️ No squeeze component found on controller');
    }
    
    // Thumbstick for zoom control instead of movement
    const thumbstickComponent = motionController.getComponent("xr-standard-thumbstick");
    if (thumbstickComponent) {
      console.log('[VR] ✅ Thumbstick component found - will use for zoom control');
      thumbstickComponent.onAxisValueChangedObservable.add((axes) => {
        // Use thumbstick Y-axis for zoom (camera distance from target)
        if (Math.abs(axes.y) > 0.1) { // Dead zone
          const camera = scene.activeCamera;
          if (camera && camera.position) {
            // Calculate zoom direction (towards/away from origin)
            const zoomDirection = camera.position.normalize();
            const zoomAmount = axes.y * 0.1; // Zoom speed
            
            // Move camera position
            camera.position.addInPlace(zoomDirection.scale(zoomAmount));
            
            // Ensure minimum distance to avoid going through the molecule
            const minDistance = 2.0;
            if (camera.position.length() < minDistance) {
              camera.position = camera.position.normalize().scale(minDistance);
            }
            
            console.log('[VR] 🔍 Zoom:', axes.y > 0 ? 'OUT' : 'IN', 'Distance:', camera.position.length().toFixed(2));
          }
        }
      });
    } else {
      console.warn('[VR] ⚠️ No thumbstick component found on controller');
    }
    
    // Set up continuous camera drag update during trigger hold
    const cameraUpdateInterval = setInterval(() => {
      if (isDragging && lastControllerPosition) {
        updateCameraFromDrag(motionController, scene);
      }
    }, 16); // ~60fps
    
    // Clean up interval when controller is removed
    motionController.onDisposeObservable.addOnce(() => {
      clearInterval(cameraUpdateInterval);
    });
    
    console.log('[VR] ✅ Controller interaction setup complete with trigger + drag camera control');
    
  } catch (error) {
    console.error('[VR] ❌ Error setting up controller interaction:', error);
  }
}

function performBondSelection(motionController, scene) {
  console.log('[VR] 🎯 performBondSelection called');
  console.log('[VR] Controller details:', {
    profile: motionController.profile,
    rootMesh: !!motionController.rootMesh,
    mesh: !!motionController.mesh,
    inputSource: !!motionController.inputSource
  });
  
  try {
    // Check if motion controller is available
    if (!motionController) {
      console.warn('[VR] ❌ Motion controller not available for bond selection');
      return;
    }
    
    // Method 1: Use WebXR input source directly
    if (motionController.inputSource && motionController.inputSource.targetRaySpace) {
      console.log('[VR] Method 1 - Using WebXR targetRaySpace');
      
      // Get the WebXR session and frame
      const xrHelper = window.vrHelper;
      if (xrHelper && xrHelper.baseExperience.sessionManager.currentFrame) {
        const frame = xrHelper.baseExperience.sessionManager.currentFrame;
        const referenceSpace = xrHelper.baseExperience.sessionManager.referenceSpace;
        
        if (frame && referenceSpace) {
          try {
            const pose = frame.getPose(motionController.inputSource.targetRaySpace, referenceSpace);
            if (pose) {
              const matrix = pose.transform.matrix;
              
              // Extract position and direction from WebXR pose matrix
              const position = new BABYLON.Vector3(matrix[12], matrix[13], matrix[14]);
              const direction = new BABYLON.Vector3(-matrix[8], -matrix[9], -matrix[10]); // Forward is -Z
              
              console.log('[VR] WebXR Controller position:', position);
              console.log('[VR] WebXR Controller direction:', direction);
              
              const ray = new BABYLON.Ray(position, direction);
              performRaycast(ray, scene);
              return;
            }
          } catch (poseError) {
            console.warn('[VR] WebXR pose extraction failed:', poseError.message);
          }
        }
      }
    }
    
    // Method 2: Try to get controller mesh for raycasting
    let controllerMesh = motionController.rootMesh;
    console.log('[VR] Method 2 - rootMesh:', !!controllerMesh);
    
    if (!controllerMesh && motionController.getComponentOfType) {
      // Method 3: Try to get pointer component
      const pointer = motionController.getComponentOfType("pointer");
      console.log('[VR] Method 3 - pointer component:', !!pointer);
      if (pointer && pointer.mesh) {
        controllerMesh = pointer.mesh;
        console.log('[VR] Method 3 - pointer.mesh:', !!controllerMesh);
      }
    }
    
    if (!controllerMesh && motionController.mesh) {
      // Method 4: Try direct mesh property
      controllerMesh = motionController.mesh;
      console.log('[VR] Method 4 - direct mesh:', !!controllerMesh);
    }
    
    // Method 5: Use controller mesh if available
    if (controllerMesh) {
      console.log('[VR] Method 5 - Creating ray from controller mesh');
      const controllerPosition = controllerMesh.absolutePosition || controllerMesh.position || new BABYLON.Vector3(0, 1.6, 0);
      
      // Get the forward direction of the controller
      let controllerDirection = new BABYLON.Vector3(0, 0, -1);
      if (controllerMesh.forward) {
        controllerDirection = controllerMesh.forward;
      } else if (controllerMesh.getDirection) {
        controllerDirection = controllerMesh.getDirection(BABYLON.Vector3.Forward());
      }
      
      console.log('[VR] Controller position:', controllerPosition);
      console.log('[VR] Controller direction:', controllerDirection);
      
      const ray = new BABYLON.Ray(controllerPosition, controllerDirection);
      performRaycast(ray, scene);
      return;
    }
    
    // Method 6: Camera fallback
    console.log('[VR] Method 6 - Using camera fallback for raycast');
    const camera = scene.activeCamera;
    if (camera) {
      console.log('[VR] Using camera position for bond selection');
      const ray = camera.getForwardRay();
      performRaycast(ray, scene);
      return;
    }
    
    console.warn('[VR] ❌ No method available for controller raycast');
    
  } catch (error) {
    console.warn('[VR] ❌ Bond selection error:', error.message);
    console.warn('[VR] Attempting camera fallback...');
    
    // Emergency fallback: use camera
    try {
      const camera = scene.activeCamera;
      if (camera) {
        const ray = camera.getForwardRay();
        performRaycast(ray, scene);
      }
    } catch (fallbackError) {
      console.error('[VR] ❌ Fallback raycast also failed:', fallbackError.message);
    }
  }
}

function performRaycast(ray, scene) {
  console.log('[VR] 🎯 Performing raycast...');
  console.log('[VR] Ray origin:', ray.origin);
  console.log('[VR] Ray direction:', ray.direction);
  
  // Add visual debug ray (optional - remove in production)
  const rayHelper = new BABYLON.RayHelper(ray);
  rayHelper.show(scene, new BABYLON.Color3(1, 1, 0)); // Yellow ray
  setTimeout(() => rayHelper.hide(), 1000); // Hide after 1 second
  
  try {
    // Count potential targets
    const bondMeshes = scene.meshes.filter(m => m.name && m.name.startsWith("bond_"));
    console.log('[VR] Found', bondMeshes.length, 'bond meshes to check');
    
    if (bondMeshes.length === 0) {
      console.warn('[VR] ⚠️ No bond meshes found in scene!');
      console.log('[VR] Available meshes:', scene.meshes.map(m => m.name).filter(Boolean));
    }
    
    // Ray cast against bond cylinders
    const pickResult = scene.pickWithRay(ray, (mesh) => {
      return mesh.name && mesh.name.startsWith("bond_");
    });
    
    console.log('[VR] Raycast result:', {
      hit: pickResult.hit,
      pickedMesh: pickResult.pickedMesh?.name,
      distance: pickResult.distance
    });
    
    if (pickResult.hit && pickResult.pickedMesh) {
      // Extract bond indices from mesh name
      const bondName = pickResult.pickedMesh.name;
      const match = bondName.match(/bond_(\d+)_(\d+)/);
      
      if (match) {
        const i = parseInt(match[1]);
        const j = parseInt(match[2]);
        
        console.log(`[VR] ✅ Selected bond ${i}-${j}`);
        
        // Visual feedback
        highlightBond(pickResult.pickedMesh);
        
        // Try to access global torsion system if available
        if (window.appState && window.appState.torsion) {
          console.log(`[VR] 🔄 Starting torsion rotation for bond ${i}-${j}`);
          window.appState.torsion.startRotation(i, j);
        } else if (window.vrTorsion) {
          console.log(`[VR] 🔄 Using VR torsion system for bond ${i}-${j}`);
          window.vrTorsion.startRotation(i, j);
        } else {
          console.warn(`[VR] ⚠️ No torsion system available - bond selected but no rotation possible`);
          console.log('[VR] Available globals:', {
            appState: !!window.appState,
            vrTorsion: !!window.vrTorsion,
            vrHelper: !!window.vrHelper
          });
        }
      } else {
        console.warn('[VR] ⚠️ Could not parse bond indices from:', bondName);
      }
    } else {
      console.log('[VR] ❌ No bond selected - ray missed all targets');
      
      // Debug: Try a general pick to see what we might hit
      const generalPick = scene.pickWithRay(ray);
      if (generalPick.hit) {
        console.log('[VR] Ray would have hit:', generalPick.pickedMesh?.name, 'at distance:', generalPick.distance);
      } else {
        console.log('[VR] Ray hit nothing at all');
      }
    }
  } catch (error) {
    console.warn('[VR] ❌ Raycast error:', error.message);
  }
}

function highlightBond(bondMesh) {
  // Create temporary highlight effect
  const originalColor = bondMesh.material.diffuseColor.clone();
  bondMesh.material.diffuseColor = new BABYLON.Color3(1, 0.5, 0); // Orange highlight
  
  setTimeout(() => {
    bondMesh.material.diffuseColor = originalColor;
  }, 500);
}

function startVRBondRotation(i, j, controller) {
  try {
    console.log(`[VR] Starting bond rotation for bond ${i}-${j}`);
    
    // Check if controller is available
    if (!controller) {
      console.warn('[VR] Controller not available for bond rotation');
      return;
    }
    
    // Use global torsion system if available
    if (window.appState && window.appState.torsion) {
      console.log(`[VR] Using global torsion system for bond ${i}-${j}`);
      window.appState.torsion.startRotation(i, j);
      
      // Set up controller tracking for rotation
      let initialControllerRotation = null;
      let isRotating = true;
      
      // Try to get controller mesh for tracking
      const controllerMesh = controller.rootMesh || controller.mesh;
      if (controllerMesh) {
        initialControllerRotation = controllerMesh.rotationQuaternion ? 
          controllerMesh.rotationQuaternion.clone() : 
          BABYLON.Quaternion.FromEulerAngles(controllerMesh.rotation.x, controllerMesh.rotation.y, controllerMesh.rotation.z);
      }
      
      // Set up rotation update mechanism
      const updateRotation = () => {
        if (!isRotating || !controllerMesh) return;
        
        const currentRotation = controllerMesh.rotationQuaternion ? 
          controllerMesh.rotationQuaternion : 
          BABYLON.Quaternion.FromEulerAngles(controllerMesh.rotation.x, controllerMesh.rotation.y, controllerMesh.rotation.z);
        
        if (initialControllerRotation) {
          // Calculate rotation delta
          const deltaRotation = currentRotation.subtract(initialControllerRotation);
          const rotationAmount = deltaRotation.y * (window.vrRotationSensitivity || 1.0);
          
          // Apply rotation to torsion system
          if (window.appState.torsion) {
            window.appState.torsion.updateRotation(rotationAmount);
          }
        }
      };
      
      // Start rotation tracking
      const rotationInterval = setInterval(updateRotation, 16); // ~60fps
      
      // Stop rotation after 5 seconds or when trigger is released
      setTimeout(() => {
        isRotating = false;
        clearInterval(rotationInterval);
        if (window.appState.torsion) {
          window.appState.torsion.stopRotation();
        }
        console.log(`[VR] Bond rotation ended for bond ${i}-${j}`);
      }, 5000);
      
    } else {
      console.warn('[VR] Global torsion system not available');
    }
    
  } catch (error) {
    console.warn('[VR] Bond rotation error:', error.message);
  }
}

function toggleVRMenu(controller, scene) {
  // Toggle floating menu panel
  console.log("[VR] Toggle menu");
  // This would show/hide the VR UI panel created in vr-ui.js
}

function updateCameraFromDrag(motionController, scene) {
  try {
    const controllerMesh = motionController.rootMesh || motionController.mesh;
    const camera = scene.activeCamera;
    
    if (!controllerMesh || !camera || !lastControllerPosition) {
      return;
    }
    
    // Get current controller position
    const currentPosition = controllerMesh.absolutePosition || controllerMesh.position;
    
    // Calculate movement delta
    const deltaPosition = currentPosition.subtract(lastControllerPosition);
    
    // Only apply movement if delta is significant enough (increased threshold)
    if (deltaPosition.length() < 0.005) {
      return;
    }
    
    // Cap maximum delta to prevent crazy jumps
    const maxDelta = 0.1;
    if (deltaPosition.length() > maxDelta) {
      console.log('[VR] 📹 Capping excessive movement delta:', deltaPosition.length().toFixed(4));
      deltaPosition.normalize();
      deltaPosition.scaleInPlace(maxDelta);
    }
    
    console.log('[VR] 📹 Camera drag update - Delta:', deltaPosition.length().toFixed(4));
    
    // Convert controller movement to camera orbital movement
    // Similar to mouse drag on desktop - horizontal movement rotates around Y axis
    // Vertical movement rotates around local X axis (up/down)
    
    const sensitivity = 0.1; // Much lower sensitivity for VR controllers
    
    // Get camera's current position relative to target (assumed to be origin)
    const cameraTarget = new BABYLON.Vector3(0, 0, 0); // Molecule center
    const cameraDirection = camera.position.subtract(cameraTarget);
    const distance = cameraDirection.length();
    
    // Convert controller delta to spherical coordinate changes
    const horizontalDelta = deltaPosition.x * sensitivity;
    const verticalDelta = deltaPosition.y * sensitivity;
    
    // Cap rotation deltas to prevent violent spinning
    const maxRotation = 0.05; // Maximum rotation per frame (~3 degrees)
    const clampedHorizontal = Math.max(-maxRotation, Math.min(maxRotation, horizontalDelta));
    const clampedVertical = Math.max(-maxRotation, Math.min(maxRotation, verticalDelta));
    
    console.log('[VR] 📹 Rotation deltas - H:', clampedHorizontal.toFixed(3), 'V:', clampedVertical.toFixed(3));
    
    // Create rotation around Y axis (horizontal movement)
    const horizontalRotation = BABYLON.Matrix.RotationY(clampedHorizontal);
    
    // Create rotation around camera's local right vector (vertical movement)
    const cameraForward = cameraDirection.normalize();
    const cameraRight = BABYLON.Vector3.Cross(cameraForward, BABYLON.Vector3.Up()).normalize();
    const verticalRotation = BABYLON.Matrix.RotationAxis(cameraRight, clampedVertical);
    
    // Apply rotations to camera direction
    let newDirection = BABYLON.Vector3.TransformCoordinates(cameraDirection, horizontalRotation);
    newDirection = BABYLON.Vector3.TransformCoordinates(newDirection, verticalRotation);
    
    // Maintain the same distance from target
    newDirection = newDirection.normalize().scale(distance);
    
    // Update camera position
    camera.position = cameraTarget.add(newDirection);
    
    // Ensure camera looks at the target
    if (camera.setTarget) {
      camera.setTarget(cameraTarget);
    }
    
    // Update last position for next frame
    lastControllerPosition = currentPosition.clone();
    
  } catch (error) {
    console.warn('[VR] Camera drag update error:', error.message);
  }
}

function handleBondSelection(bondMesh) {
  console.log('[VR] 🎯 Bond selected via pointer:', bondMesh.name);
  
  try {
    // Extract bond indices from mesh name
    const bondName = bondMesh.name;
    const match = bondName.match(/bond_(\d+)_(\d+)/);
    
    if (match) {
      const i = parseInt(match[1]);
      const j = parseInt(match[2]);
      
      console.log(`[VR] ✅ Parsed bond indices: ${i}-${j}`);
      
      // Visual feedback
      highlightBond(bondMesh);
      
      // Try to access global torsion system
      if (window.appState && window.appState.torsion) {
        console.log(`[VR] 🔄 Starting torsion rotation for bond ${i}-${j}`);
        window.appState.torsion.startRotation(i, j);
        
        // Show feedback to user
        console.log('[VR] 🎮 Bond rotation started - you can now rotate the molecule');
        
        // Optional: Add visual indicator that rotation is active
        if (bondMesh.material) {
          const originalEmissive = bondMesh.material.emissiveColor?.clone() || new BABYLON.Color3(0, 0, 0);
          bondMesh.material.emissiveColor = new BABYLON.Color3(0, 1, 0); // Green for active
          
          setTimeout(() => {
            bondMesh.material.emissiveColor = originalEmissive;
          }, 2000);
        }
        
      } else if (window.vrTorsion) {
        console.log(`[VR] 🔄 Using VR torsion system for bond ${i}-${j}`);
        window.vrTorsion.startRotation(i, j);
      } else {
        console.warn(`[VR] ⚠️ No torsion system available for rotation`);
        console.log('[VR] Available globals:', {
          appState: !!window.appState,
          vrTorsion: !!window.vrTorsion,
          vrHelper: !!window.vrHelper
        });
      }
    } else {
      console.warn('[VR] ⚠️ Could not parse bond indices from mesh name:', bondName);
    }
  } catch (error) {
    console.error('[VR] ❌ Bond selection error:', error.message);
  }
}

// Export the main setup function
export { setupVRFeatures };
