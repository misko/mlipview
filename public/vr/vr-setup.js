// vr/vr-setup.js - VR initialization and controller setup
export async function setupVR(engine, scene) {
  // Check VR support
  const xrHelper = await scene.createDefaultXRExperienceAsync({
    floorMeshes: [], // No floor mesh needed for molecular viewer
    optionalFeatures: true
  });
  
  if (!xrHelper.baseExperience) {
    console.warn("WebXR not supported");
    return null;
  }
  
  // VR-specific settings
  xrHelper.baseExperience.sessionManager.onXRSessionInit.add(() => {
    console.log("VR Session started");
    
    // Adjust lighting for VR
    scene.lights.forEach(light => {
      if (light instanceof BABYLON.HemisphericLight) {
        light.intensity = 0.8; // Slightly dimmer for VR comfort
      }
    });
    
    // Set comfortable viewing distance for molecules
    if (scene.activeCamera) {
      scene.activeCamera.setTarget(BABYLON.Vector3.Zero());
    }
  });
  
  // Hand controller setup
  xrHelper.input.onControllerAddedObservable.add((controller) => {
    console.log("VR Controller added:", controller.inputSource.handedness);
    
    controller.onMotionControllerInitObservable.add((motionController) => {
      setupControllerInteraction(motionController, scene);
    });
  });
  
  return xrHelper;
}

function setupControllerInteraction(motionController, scene) {
  // Get the main trigger component
  const triggerComponent = motionController.getComponent("xr-standard-trigger");
  const squeezeComponent = motionController.getComponent("xr-standard-squeeze");
  
  if (triggerComponent) {
    // Trigger for bond selection/rotation
    triggerComponent.onButtonStateChangedObservable.add((component) => {
      if (component.pressed) {
        // Ray casting for bond selection
        performBondSelection(motionController, scene);
      }
    });
  }
  
  if (squeezeComponent) {
    // Squeeze for menu activation
    squeezeComponent.onButtonStateChangedObservable.add((component) => {
      if (component.pressed) {
        toggleVRMenu(motionController, scene);
      }
    });
  }
  
  // Thumbstick for rotation sensitivity
  const thumbstickComponent = motionController.getComponent("xr-standard-thumbstick");
  if (thumbstickComponent) {
    thumbstickComponent.onAxisValueChangedObservable.add((axes) => {
      // Use thumbstick Y-axis to control rotation sensitivity
      const sensitivity = 1 + axes.y; // 0.5x to 2x speed
      window.vrRotationSensitivity = Math.max(0.1, sensitivity);
    });
  }
}

function performBondSelection(motionController, scene) {
  // Get controller ray
  const ray = new BABYLON.Ray(
    motionController.pointer.position,
    motionController.pointer.forward
  );
  
  // Ray cast against bond cylinders
  const pickResult = scene.pickWithRay(ray, (mesh) => {
    return mesh.name && mesh.name.startsWith("bond_");
  });
  
  if (pickResult.hit && pickResult.pickedMesh) {
    // Extract bond indices from mesh name
    const bondName = pickResult.pickedMesh.name;
    const match = bondName.match(/bond_(\d+)_(\d+)/);
    
    if (match) {
      const i = parseInt(match[1]);
      const j = parseInt(match[2]);
      
      // Trigger bond rotation (this would connect to existing torsion system)
      console.log(`VR: Selected bond ${i}-${j}`);
      
      // Visual feedback
      highlightBond(pickResult.pickedMesh);
      
      // Start rotation mode
      startVRBondRotation(i, j, motionController);
    }
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
  // This would integrate with the existing torsion rotation system
  // Store initial controller position/rotation for calculating rotation delta
  
  const initialRotation = controller.pointer.rotationQuaternion.clone();
  let isRotating = true;
  
  const rotationObserver = controller.onMotionControllerInitObservable.add(() => {
    if (!isRotating) return;
    
    // Calculate rotation delta
    const currentRotation = controller.pointer.rotationQuaternion;
    const deltaRotation = currentRotation.subtract(initialRotation);
    
    // Convert to angle and apply to bond
    const angleDelta = deltaRotation.toEulerAngles().y * 180 / Math.PI;
    const sensitivity = window.vrRotationSensitivity || 1;
    
    // This would call the existing rotation system
    // applyBondRotation(i, j, angleDelta * sensitivity);
  });
  
  // Stop rotation on trigger release
  const triggerComponent = controller.motionController.getComponent("xr-standard-trigger");
  triggerComponent.onButtonStateChangedObservable.addOnce((component) => {
    if (!component.pressed) {
      isRotating = false;
      controller.onMotionControllerInitObservable.remove(rotationObserver);
    }
  });
}

function toggleVRMenu(controller, scene) {
  // Toggle floating menu panel
  console.log("VR: Toggle menu");
  // This would show/hide the VR UI panel created in vr-ui.js
}
