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
  
  // Lite mode: only expose energy value (no control buttons)
  return {
    advancedTexture,
    energyValue
  };
}
