// vr/vr-ui.js - VR-compatible UI using Babylon.js GUI
export function createVRUI(scene) {
  // Create a GUI texture for VR
  const advancedTexture = BABYLON.GUI.AdvancedDynamicTexture.CreateFullscreenUI("VR_UI", false, scene);
  // Improve readability/stability in XR: render at a consistent ideal size
  try {
      advancedTexture.idealHeight = 1024; // crisp text in HMD
      advancedTexture.renderAtIdealSize = true;
      advancedTexture.useInvalidateRectOptimization = false; // avoid partial redraw glitches in XR
  } catch {}
  
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
  bondBar.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_BOTTOM;
  bondBar.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  bondBar.width = "90%";
  bondBar.isVisible = true; // always show; label will instruct when nothing selected
  bondBar.paddingBottom = "24px"; // keep safely inside view
  try { bondBar.zIndex = 1000; } catch {}
  advancedTexture.addControl(bondBar);
  try { energyPanel.zIndex = 500; } catch {}
  advancedTexture.addControl(energyPanel);

  const row = new BABYLON.GUI.StackPanel();
  row.isVertical = false;
  row.width = "100%";
  row.height = "100%";
  row.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  row.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_CENTER;
  row.spacing = 10;
  bondBar.addControl(row);

  function mkBtn(text, w = "120px") {
    const b = BABYLON.GUI.Button.CreateSimpleButton("", text);
    b.height = "100%";
    b.width = w;
    b.color = "#e9f6ff";
    b.fontSize = 30;
    b.thickness = 0;
    b.paddingLeft = "8px";
    b.paddingRight = "8px";
    b.background = "rgba(35,45,60,0.9)";
    b.cornerRadius = 10;
    b.onPointerEnterObservable.add(() => b.background = "rgba(50,65,85,0.95)");
    b.onPointerOutObservable.add(() => b.background = "rgba(35,45,60,0.9)");
    return b;
  }

  const btnMinus = mkBtn("⟲ −");
  const labelRect = new BABYLON.GUI.Rectangle("bondLblRect");
  labelRect.height = "100%";
  labelRect.width = "260px";
  labelRect.thickness = 0;
  labelRect.background = "rgba(25,32,44,0.9)";
  labelRect.cornerRadius = 10;
  const bondLabel = new BABYLON.GUI.TextBlock("bondLbl", "Select a bond to rotate");
  bondLabel.color = "#d7e6ff";
  bondLabel.fontSize = 30;
  bondLabel.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  bondLabel.textVerticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_CENTER;
  labelRect.addControl(bondLabel);
  const btnPlus  = mkBtn("⟲ +");
  const btnSide  = mkBtn("Side: j");
  const btnRec   = mkBtn("recompute: off");

  row.addControl(btnMinus);
  row.addControl(labelRect);
  row.addControl(btnPlus);
  row.addControl(btnSide);
  row.addControl(btnRec);

  // Note: placing btnRec overlapping step+ if needed; adjust columns later if you want separate columns

  const bondUIState = { side: 'j', step: 5, recompute: false };

  btnSide.onPointerUpObservable.add(() => {
    bondUIState.side = bondUIState.side === 'j' ? 'i' : 'j';
    btnSide.textBlock.text = `Side: ${bondUIState.side}`;
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
      btnRec,
      state: bondUIState
    }
  };
}

// Create a camera-anchored HUD for XR sessions using a mesh in front of the camera
export function createVRUIOnCamera(scene, xrCamera) {
  // Runtime overrides for quick tuning on device
  const cfg = (typeof window !== 'undefined' && window.vrHudConfig) ? window.vrHudConfig : {};
  const fontScale = cfg.fontScale || 1.2;
  const pos = cfg.position || { x: 0, y: -0.15, z: 1.1 };

  // Root transform, parented to camera
  const hudRoot = new BABYLON.TransformNode('vrHudRoot', scene);
  hudRoot.parent = xrCamera;
  hudRoot.position = new BABYLON.Vector3(pos.x, pos.y, pos.z);

  // Utility to make a button plane with its own ADT (only area pickable)
  function mkButtonPlane(name, w, h, localPos, text, fontPx) {
    const mesh = BABYLON.MeshBuilder.CreatePlane(name, { width: w, height: h, sideOrientation: BABYLON.Mesh.DOUBLESIDE }, scene);
    mesh.isPickable = true; // only this area blocks picks
    mesh.renderingGroupId = 3;
    mesh.parent = hudRoot;
    mesh.position = localPos.clone();
    const adt = BABYLON.GUI.AdvancedDynamicTexture.CreateForMesh(mesh, Math.round(w * 1800), Math.round(h * 1800), false);
    try { adt.useInvalidateRectOptimization = false; } catch {}
    const btn = BABYLON.GUI.Button.CreateSimpleButton(name + '_btn', text);
    btn.width = 1;
    btn.height = 1;
    btn.color = '#e9f6ff';
    btn.fontSize = Math.round((fontPx || 34) * fontScale);
    btn.thickness = 0;
    btn.background = 'rgba(35,45,60,0.9)';
    btn.cornerRadius = 10;
    btn.onPointerEnterObservable.add(() => btn.background = 'rgba(50,65,85,0.95)');
    btn.onPointerOutObservable.add(() => btn.background = 'rgba(35,45,60,0.9)');
    adt.addControl(btn);
    return { mesh, adt, btn };
  }

  // Non-interactive label plane (pass-through)
  function mkLabelPlane(name, w, h, localPos, initialText, fontPx) {
    const mesh = BABYLON.MeshBuilder.CreatePlane(name, { width: w, height: h, sideOrientation: BABYLON.Mesh.DOUBLESIDE }, scene);
    mesh.isPickable = false; // never block picks
    mesh.renderingGroupId = 3;
    mesh.parent = hudRoot;
    mesh.position = localPos.clone();
    const adt = BABYLON.GUI.AdvancedDynamicTexture.CreateForMesh(mesh, Math.round(w * 1800), Math.round(h * 1800), false);
    try { adt.useInvalidateRectOptimization = false; } catch {}
    const rect = new BABYLON.GUI.Rectangle(name + '_rect');
    rect.thickness = 0;
    rect.background = 'rgba(25,32,44,0.9)';
    rect.cornerRadius = 10;
    adt.addControl(rect);
    const txt = new BABYLON.GUI.TextBlock(name + '_text', initialText);
    txt.color = '#d7e6ff';
    txt.fontSize = Math.round((fontPx || 34) * fontScale);
    txt.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
    txt.textVerticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_CENTER;
    rect.addControl(txt);
    return { mesh, adt, rect, txt };
  }

  // Small energy panel (non-pickable) at top-right of HUD root
  const energyW = 0.32, energyH = 0.14;
  const energy = mkLabelPlane('energyPanelXR', energyW, energyH, new BABYLON.Vector3(0.35, 0.22, 0), '0.000', 24);
  // Add a title above energy value
  const energyTitle = new BABYLON.GUI.TextBlock('energyTitleXR', 'Mock Energy');
  energyTitle.color = '#a4b0c0';
  energyTitle.fontSize = 18;
  energyTitle.top = '-40px';
  try { energy.rect.addControl(energyTitle); } catch {}
  const energyValue = energy.txt;

  // Layout sizes in meters
  const hBtn = 0.10; // 10 cm
  const wBtn = 0.22; // 22 cm
  const wLabel = 0.62;
  const wRec = 0.34;
  const spacing = 0.03;
  const rowDY = 0.14; // vertical distance between rows
  const row1Y = 0.06;
  const row2Y = -0.08;

  // Compute positions
  const minusX = -(wLabel / 2 + spacing + wBtn / 2);
  const plusX  = +(wLabel / 2 + spacing + wBtn / 2);
  const sideRecTotal = wBtn + spacing + wRec;
  const startX = -sideRecTotal / 2 + wBtn / 2;
  const sideX = startX;
  const recX = startX + wBtn / 2 + spacing + wRec / 2;

  const minus = mkButtonPlane('hudMinus', wBtn, hBtn, new BABYLON.Vector3(minusX, row1Y, 0), '⟲ −');
  const label = mkLabelPlane('hudLabel', wLabel, hBtn, new BABYLON.Vector3(0, row1Y, 0), 'Select a bond to rotate');
  const plus  = mkButtonPlane('hudPlus',  wBtn, hBtn, new BABYLON.Vector3(plusX, row1Y, 0), '⟲ +');
  const side  = mkButtonPlane('hudSide',  wBtn, hBtn, new BABYLON.Vector3(sideX, row2Y, 0), 'Side: j');
  const rec   = mkButtonPlane('hudRecompute', wRec, hBtn, new BABYLON.Vector3(recX, row2Y, 0), 'recompute: off');

  // State shared with app
  const bondUIState = { side: 'j', step: 5, recompute: false };

  // GUI press guard (only when real button planes pressed)
  const setGuiActive = (v) => {
    try {
      window.vrGuiPointerActive = !!v;
      if (v) {
        if (console && console.debug) console.debug('[HUD] GUI active: down');
        window.vrGuiActiveUntil = 0;
      } else {
        const until = Date.now() + 150;
        window.vrGuiActiveUntil = until;
        if (console && console.debug) console.debug('[HUD] GUI inactive: up, grace until', new Date(until).toLocaleTimeString());
      }
    } catch {}
  };
  const wirePressGuards = (btn) => {
    if (!btn) return;
    btn.onPointerDownObservable.add(() => setGuiActive(true));
    btn.onPointerUpObservable.add(() => setTimeout(() => setGuiActive(false), 0));
  };
  [minus.btn, plus.btn, side.btn, rec.btn].forEach(wirePressGuards);

  // Side/recompute handlers
  side.btn.onPointerUpObservable.add(() => {
    bondUIState.side = bondUIState.side === 'j' ? 'i' : 'j';
    side.btn.textBlock.text = `Side: ${bondUIState.side}`;
  });
  rec.btn.onPointerUpObservable.add(() => {
    bondUIState.recompute = !bondUIState.recompute;
    rec.btn.textBlock.text = `recompute: ${bondUIState.recompute ? 'on' : 'off'}`;
  });

  return {
    advancedTexture: null, // per-mesh ADTs; no single texture
    rootMesh: hudRoot,
    energyValue,
    bond: {
      label: label.txt,
      btnMinus: minus.btn,
      btnPlus: plus.btn,
      btnSide: side.btn,
      btnRec: rec.btn,
      state: bondUIState
    }
  };
}
