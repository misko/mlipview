// vr/vr-ui.js - VR-compatible UI using Babylon.js GUI
export function createVRUI(scene) {
  // Create a GUI texture for VR
  const advancedTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI("VR_UI", true, scene);
  
  // Energy display panel
  const energyPanel = new BABYLON.GUI.Rectangle("energyPanel");
  energyPanel.width = "300px";
  energyPanel.height = "120px";
  energyPanel.cornerRadius = 10;
  energyPanel.color = "rgba(255,255,255,0.08)";
  energyPanel.background = "rgba(15,18,24,0.9)";
  energyPanel.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_RIGHT;
  energyPanel.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_TOP;
  energyPanel.top = "20px";
  energyPanel.left = "-20px";
  advancedTexture.addControl(energyPanel);
  
  // Energy value text
  const energyTitle = new BABYLON.GUI.TextBlock("energyTitle", "Mock Energy");
  energyTitle.color = "#a4b0c0";
  energyTitle.fontSize = "14px";
  energyTitle.top = "-30px";
  energyPanel.addControl(energyTitle);
  
  const energyValue = new BABYLON.GUI.TextBlock("energyValue", "0.000");
  energyValue.color = "#f5ffd1";
  energyValue.fontSize = "24px";
  energyValue.fontWeight = "bold";
  energyValue.top = "5px";
  energyPanel.addControl(energyValue);
  
  // Control buttons panel
  const controlPanel = new BABYLON.GUI.StackPanel("controlPanel");
  controlPanel.width = "300px";
  controlPanel.height = "200px";
  controlPanel.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_RIGHT;
  controlPanel.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_TOP;
  controlPanel.top = "160px";
  controlPanel.left = "-20px";
  advancedTexture.addControl(controlPanel);
  
  // Create VR-friendly buttons
  const createVRButton = (text, callback) => {
    const button = BABYLON.GUI.Button.CreateSimpleButton(`btn_${text}`, text);
    button.width = "280px";
    button.height = "40px";
    button.color = "#f5ffd1";
    button.background = "rgba(255,255,255,0.1)";
    button.cornerRadius = 5;
    button.fontSize = "16px";
    button.paddingTop = "5px";
    button.paddingBottom = "5px";
    
    button.onPointerClickObservable.add(callback);
    
    // VR hover effects
    button.onPointerEnterObservable.add(() => {
      button.background = "rgba(255,255,255,0.2)";
    });
    button.onPointerOutObservable.add(() => {
      button.background = "rgba(255,255,255,0.1)";
    });
    
    return button;
  };
  
  const forcesBtn = createVRButton("Forces: ON", () => {});
  const lengthBtn = createVRButton("Length: Normalized", () => {});
  const plotBtn = createVRButton("Plot: ON", () => {});
  const rebuildBtn = createVRButton("Rebuild from rotations", () => {});
  const exportBtn = createVRButton("Export XYZ", () => {});
  
  controlPanel.addControl(forcesBtn);
  controlPanel.addControl(lengthBtn);
  controlPanel.addControl(plotBtn);
  controlPanel.addControl(rebuildBtn);
  controlPanel.addControl(exportBtn);
  
  return {
    advancedTexture,
    energyValue,
    buttons: {
      forces: forcesBtn,
      length: lengthBtn,
      plot: plotBtn,
      rebuild: rebuildBtn,
      export: exportBtn
    }
  };
}
