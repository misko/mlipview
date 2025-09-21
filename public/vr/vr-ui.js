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
  // Bond control bottom bar (2D overlay)
  const bondBar = new BABYLON.GUI.Rectangle("bondBar");
  bondBar.height = "80px";
  bondBar.thickness = 0;
  bondBar.background = "rgba(15,18,24,0.66)";
  bondBar.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_BOTTOM;
  bondBar.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  bondBar.width = 0.9;
  bondBar.isVisible = false; // show only when a bond is selected
  advancedTexture.addControl(bondBar);

  const bondGrid = new BABYLON.GUI.Grid();
  // 8 columns: [minus, label, plus, side, step-, step, step+, recompute]
  for (let i = 0; i < 8; i++) bondGrid.addColumnDefinition(1/8);
  bondGrid.addRowDefinition(1.0);
  bondBar.addControl(bondGrid);

  function mkBtn(text, w = "auto") {
    const b = BABYLON.GUI.Button.CreateSimpleButton("", text);
    b.height = 1.0;
    b.width = w;
    b.color = "#e9f6ff";
    b.fontSize = 30;
    b.thickness = 0;
    b.paddingLeft = "6px";
    b.paddingRight = "6px";
    b.background = "rgba(35,45,60,0.9)";
    b.cornerRadius = 10;
    b.onPointerEnterObservable.add(() => b.background = "rgba(50,65,85,0.95)");
    b.onPointerOutObservable.add(() => b.background = "rgba(35,45,60,0.9)");
    return b;
  }

  const btnMinus = mkBtn("⟲ −");
  const bondLabel = new BABYLON.GUI.TextBlock("bondLbl", "Bond");
  bondLabel.color = "#d7e6ff";
  bondLabel.fontSize = 30;
  bondLabel.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  const btnPlus  = mkBtn("⟲ +");
  const btnSide  = mkBtn("Side: j");
  const btnStepM = mkBtn("Step −");
  const stepText = new BABYLON.GUI.TextBlock("stepTxt", "Step: 5°");
  stepText.color = "#d7e6ff";
  stepText.fontSize = 28;
  stepText.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  const btnStepP = mkBtn("Step +");
  const btnRec   = mkBtn("recompute: off");

  bondGrid.addControl(btnMinus, 0, 0);
  bondGrid.addControl(bondLabel, 0, 1);
  bondGrid.addControl(btnPlus,  0, 2);
  bondGrid.addControl(btnSide,  0, 3);
  bondGrid.addControl(btnStepM, 0, 4);
  bondGrid.addControl(stepText, 0, 5);
  bondGrid.addControl(btnStepP, 0, 6);
  bondGrid.addControl(btnRec,   0, 7);

  // Note: placing btnRec overlapping step+ if needed; adjust columns later if you want separate columns

  const bondUIState = { side: 'j', step: 5, recompute: false };

  btnSide.onPointerUpObservable.add(() => {
    bondUIState.side = bondUIState.side === 'j' ? 'i' : 'j';
    btnSide.textBlock.text = `Side: ${bondUIState.side}`;
  });
  btnStepM.onPointerUpObservable.add(() => {
    bondUIState.step = Math.max(1, Math.min(45, bondUIState.step - 1));
    stepText.text = `Step: ${bondUIState.step}°`;
  });
  btnStepP.onPointerUpObservable.add(() => {
    bondUIState.step = Math.max(1, Math.min(45, bondUIState.step + 1));
    stepText.text = `Step: ${bondUIState.step}°`;
  });
  btnRec.onPointerUpObservable.add(() => {
    bondUIState.recompute = !bondUIState.recompute;
    btnRec.textBlock.text = `recompute: ${bondUIState.recompute ? 'on' : 'off'}`;
  });

  return {
    advancedTexture,
    energyValue,
    // overlay bond UI
    bond: {
      bar: bondBar,
      label: bondLabel,
      btnMinus,
      btnPlus,
      btnSide,
      btnStepM,
      btnStepP,
      btnRec,
      stepText,
      state: bondUIState
    }
  };
}
