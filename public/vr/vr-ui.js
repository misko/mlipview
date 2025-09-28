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
  // Move overlay energy panel lower by an additional ~10% of screen height
  energyPanel.top = "30%";
  energyPanel.left = "-20px";
  advancedTexture.addControl(energyPanel);
  
  // Energy value text
  const energyTitle = new BABYLON.GUI.TextBlock("energyTitle", "Mock Energy");
  energyTitle.color = "#a4b0c0";
  energyTitle.fontSize = "18px";
  energyTitle.fontWeight = "bold";
  energyTitle.top = "-30px";
  energyPanel.addControl(energyTitle);
  
  const energyValue = new BABYLON.GUI.TextBlock("energyValue", "0.000");
  energyValue.color = "#f5ffd1";
  energyValue.fontSize = "30px";
  energyValue.fontWeight = "bold";
  energyValue.top = "5px";
  energyPanel.addControl(energyValue);
  
  // Lite mode: only expose energy value (no control buttons)
  // Bond control bottom bar (2D overlay)
  const bondBar = new BABYLON.GUI.Rectangle("bondBar");
  bondBar.height = "240px";
  bondBar.thickness = 0;
  bondBar.background = "rgba(15,18,24,0.66)";
  bondBar.verticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_BOTTOM;
  bondBar.horizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
  bondBar.width = "90%";
  bondBar.isVisible = true; // always show; label will instruct when nothing selected
  bondBar.paddingBottom = "24px"; // keep safely inside view
  try { bondBar.zIndex = 1000; } catch {}
  advancedTexture.addControl(bondBar);
  try { energyPanel.zIndex = 500; } catch {}

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
    // Triple overlay button font size for readability
    b.fontSize = 108;
    b.fontWeight = "bold";
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
  const btnPlus  = mkBtn("⟲ +");
  // Side toggle removed; rotation direction now determined by bond orientation cycle (selection logic)
  row.addControl(btnPlus);
  row.addControl(btnMinus);
  const bondUIState = { step: 5 };
  // recompute toggle removed

  return {
    advancedTexture,
    energyValue,
    // overlay bond UI
    bond: {
      bar: bondBar,
  btnMinus,
  btnPlus,
  state: bondUIState
    },
    dispose() {
      try { advancedTexture.dispose(); } catch {}
    }
  };
}

// Create a camera-anchored HUD for XR sessions using a mesh in front of the camera
export function createVRUIOnCamera(scene, xrCamera) {
  // Runtime overrides for quick tuning on device
  const cfg = (typeof window !== 'undefined' && window.vrHudConfig) ? window.vrHudConfig : {};
  const fontScale = cfg.fontScale || 1.2;
  // Position the HUD, shifted lower by an additional ~10% for comfort
  // Lower the entire HUD noticeably; use runtime override if provided
  const pos = cfg.position || { x: 0, y: -0.06, z: 1.1 };

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
    const fs = Math.round((fontPx || 42) * fontScale);
    // Set font size on both the Button and its inner TextBlock to ensure it takes effect
    btn.fontSize = fs;
    try {
      if (btn.textBlock) {
        btn.textBlock.fontWeight = 'bold';
        btn.textBlock.fontSize = fs;
      }
    } catch {}
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
    txt.fontSize = Math.round((fontPx || 38) * fontScale);
  try { txt.fontWeight = 'bold'; } catch {}
    txt.textHorizontalAlignment = BABYLON.GUI.Control.HORIZONTAL_ALIGNMENT_CENTER;
    txt.textVerticalAlignment = BABYLON.GUI.Control.VERTICAL_ALIGNMENT_CENTER;
    rect.addControl(txt);
    return { mesh, adt, rect, txt };
  }

  // Small energy panel (non-pickable) at top-right of HUD root (near top of FOV)
  const energyW = 0.40, energyH = 0.16;
  // Energy panel local position (root shift handles 10% lowering globally)
  // Swap energy to the left side and raise higher in the HUD
  const energyPos = (cfg.energyPos && typeof cfg.energyPos === 'object')
    ? new BABYLON.Vector3(cfg.energyPos.x ?? -0.45, cfg.energyPos.y ?? 0.26, cfg.energyPos.z ?? 0)
    : new BABYLON.Vector3(-0.45, 0.26, 0);
  const energy = mkLabelPlane('energyPanelXR', energyW, energyH, energyPos, '0.000', 45);
  // Add a title above energy value
  const energyTitle = new BABYLON.GUI.TextBlock('energyTitleXR', 'Mock Energy');
  energyTitle.color = '#a4b0c0';
  energyTitle.fontSize = 33;
  energyTitle.fontWeight = 'bold';
  energyTitle.top = '-40px';
  try { energy.rect.addControl(energyTitle); } catch {}
  const energyValue = energy.txt;
  try { energyValue.fontWeight = 'bold'; } catch {}

  // ---- Mini energy vs step plot (top-left) ----
  function createHudPlot(scene, parent, opts = {}) {
    const w = opts.width || 0.60;   // meters
    const h = opts.height || 0.24;  // meters
    // Plot local position; lower it and keep the 10% right shift. Allow override via cfg.plotPos
    // Swap plot to the right side and raise higher in the HUD
    const defaultPlotPos = new BABYLON.Vector3(0.39, 0.30, 0.001);
    const origin = opts.position || (cfg.plotPos
      ? new BABYLON.Vector3(cfg.plotPos.x ?? 0.39, cfg.plotPos.y ?? 0.30, cfg.plotPos.z ?? 0.001)
      : defaultPlotPos);
    const maxPoints = opts.maxPoints || 100;
    const root = new BABYLON.TransformNode('hudPlotRoot', scene);
    root.parent = parent;
    root.position = origin.clone();

    // Background plane with dynamic texture for plot
  // Create a plane rendering only its front face and make it always face the camera
  const bg = BABYLON.MeshBuilder.CreatePlane('hudPlotBg', { width: w, height: h, sideOrientation: BABYLON.Mesh.FRONTSIDE }, scene);
    bg.parent = root;
    bg.position = new BABYLON.Vector3(0, 0, 0);
  // Make the plane always face the camera to avoid mirrored/back-side views
  try { bg.rotation = new BABYLON.Vector3(0, 0, 0); } catch {}
  try { bg.billboardMode = BABYLON.AbstractMesh.BILLBOARDMODE_ALL; } catch {}
    bg.isPickable = false;
    bg.renderingGroupId = 3;
  const mat = new BABYLON.StandardMaterial('hudPlotMat', scene);
  mat.disableLighting = true;
  // Cull back face so we never see mirrored text
  mat.backFaceCulling = true;
    mat.alpha = 0.98;
    // Draw onto a high-res dynamic texture for crisp lines
    const texW = 1024, texH = 512;
  const dyn = new BABYLON.DynamicTexture('hudPlotDT', { width: texW, height: texH }, scene, false);
    dyn.hasAlpha = true;
    // Fix orientation: invert both axes so the canvas draws appear non-mirrored in world
    try {
      // Keep vertical flip so bottom of canvas maps to bottom in world
      dyn.vScale = -1; dyn.vOffset = 1;
      // Do not flip horizontally so left stays left (Y-axis/"Energy" on left)
      dyn.uScale = 1; dyn.uOffset = 0;
    } catch {}
    mat.emissiveTexture = dyn;
    mat.diffuseTexture = dyn;
    bg.material = mat;

  // Data buffers: energies and corresponding step indices
  // We decouple the X axis (step index) from wall-clock time so the plot
  // reflects user-driven structural changes / sampling steps, matching
  // desktop semantics.
  let data = [];
  let steps = [];
  let minY = Infinity, maxY = -Infinity;

    function redraw() {
      const ctx = dyn.getContext();
      // Clear
      ctx.clearRect(0, 0, texW, texH);
      // Draw in normal orientation (no canvas flipping)
      // Background
      ctx.fillStyle = 'rgba(26, 32, 44, 0.92)';
      ctx.fillRect(0, 0, texW, texH);
      // Axes
      ctx.strokeStyle = 'rgba(200, 200, 200, 0.7)';
      ctx.lineWidth = 2;
      // Padding inside texture
      const padL = 60, padR = 20, padT = 20, padB = 40;
      const x0 = padL, x1 = texW - padR, y0 = texH - padB, y1 = padT;
      // X-axis
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x1, y0);
      ctx.stroke();
      // Y-axis
      ctx.beginPath();
      ctx.moveTo(x0, y0);
      ctx.lineTo(x0, y1);
      ctx.stroke();
      // Axis labels (bold)
    ctx.fillStyle = 'rgba(220, 230, 255, 0.95)';
  ctx.font = 'bold 32px sans-serif';
  // X label: Step (bottom center) - with flipped texture, this maps to bottom in world
  const xLabel = 'Step';
  const xLabelWidth = ctx.measureText(xLabel).width;
  ctx.fillText(xLabel, x0 + (x1 - x0) / 2 - xLabelWidth / 2, texH - 14);
      // Y label: Energy (rotated)
      ctx.save();
  // Place the Y label slightly to the left of the axis
  ctx.translate(14, y0 - (y0 - y1) / 2);
      ctx.rotate(-Math.PI / 2);
      const yLabel = 'Energy';
      const yLabelWidth = ctx.measureText(yLabel).width;
      ctx.fillText(yLabel, -yLabelWidth / 2, 0);
      ctx.restore();
      // Plot line
      if (data.length >= 2) {
        minY = Math.min(...data);
        maxY = Math.max(...data);
        if (minY === maxY) { minY -= 1; maxY += 1; }
        const n = data.length;
        // Determine step range for scaling X properly (even if steps not contiguous)
        const minStep = steps[0];
        const maxStep = steps[steps.length - 1];
        const span = Math.max(1, maxStep - minStep);
        ctx.strokeStyle = 'rgba(245, 255, 209, 1)';
        ctx.lineWidth = 4;
        ctx.beginPath();
        for (let i = 0; i < n; i++) {
          const sNorm = (steps[i] - minStep) / span; // 0..1
          const x = x0 + sNorm * (x1 - x0);
          const ny = (data[i] - minY) / (maxY - minY);
          const y = y0 - ny * (y0 - y1);
          if (i === 0) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        }
        ctx.stroke();
      }
      dyn.update(false);
    }

    function addPoint(v, stepIndex) {
      // Accept only finite energies
      if (!Number.isFinite(v)) return;
      const nextStep = (typeof stepIndex === 'number' ? stepIndex : (steps.length ? steps[steps.length - 1] + 1 : 0));
      data.push(v);
      steps.push(nextStep);
      if (data.length > maxPoints) { data.shift(); steps.shift(); }
      redraw();
    }
    function reset() { data = []; steps = []; redraw(); }
    function setVisible(visible) { bg.setEnabled(!!visible); }

    // Initial draw
    redraw();

  return { addPoint, reset, setVisible, root, mesh: bg, mat, dyn };
  }

  const hudPlot = createHudPlot(scene, hudRoot);

  // Layout sizes in meters (reverted physical size; we scale font instead)
  const hBtn = 0.10; // 10 cm
  const wBtn = 0.22; // 22 cm
  // recompute removed; no wide button
  const spacing = 0.03;
  const rowDY = 0.16; // vertical distance between rows
  // Place buttons lower at bottom area of HUD (below center), single row now
  // Lower further by default; allow runtime override via window.vrHudConfig.buttonRowY
  const rowY = (typeof cfg.buttonRowY === 'number') ? cfg.buttonRowY : -0.32;

  // Compute positions
  // Single row positions now: +  -  (center formerly side toggle)
  const plusX  = - (wBtn + spacing*0.5);
  const minusX = + (wBtn + spacing*0.5);
  const molsX  = 0; // Will place molecules button in a second row below center
  // recX removed

  // Triple font size for XR control buttons only via local fontPx
  const minus = mkButtonPlane('hudMinus', wBtn, hBtn, new BABYLON.Vector3(minusX, rowY, 0), '⟲ −', 42 * 3);
  const plus  = mkButtonPlane('hudPlus',  wBtn, hBtn, new BABYLON.Vector3(plusX, rowY, 0), '⟲ +', 42 * 3);
  const mols  = mkButtonPlane('hudMols',  wBtn*1.4, hBtn, new BABYLON.Vector3(molsX, rowY - (rowDY*1.1), 0), 'Molecules', 42 * 2.2);
  // recompute plane removed

  // State shared with app
  const bondUIState = { step: 5 }; // side removed; orientation from selection model

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
  [minus.btn, plus.btn, mols.btn].forEach(wirePressGuards);

  // Side toggle removed (orientation handled via selection clicks)

  return {
    advancedTexture: null, // per-mesh ADTs; no single texture
    rootMesh: hudRoot,
    energyValue,
    plot: hudPlot,
    bond: {
      btnMinus: minus.btn,
      btnPlus: plus.btn,
      btnSide: null,
      state: bondUIState
    },
    btnMolecules: mols.btn,
    dispose() {
      // Dispose GUI textures first
      try { minus.adt && minus.adt.dispose(); } catch {}
      try { plus.adt && plus.adt.dispose(); } catch {}
      try { mols.adt && mols.adt.dispose(); } catch {}
      try { energy.adt && energy.adt.dispose(); } catch {}
      // Dispose plot textures/materials/mesh
      try { hudPlot && hudPlot.dyn && hudPlot.dyn.dispose(); } catch {}
      try { hudPlot && hudPlot.mat && hudPlot.mat.dispose(); } catch {}
      try { hudPlot && hudPlot.mesh && hudPlot.mesh.dispose(); } catch {}
      // Dispose button/label meshes
      try { minus.mesh && minus.mesh.dispose(); } catch {}
      try { plus.mesh && plus.mesh.dispose(); } catch {}
      try { mols.mesh && mols.mesh.dispose(); } catch {}
      try { energy.mesh && energy.mesh.dispose(); } catch {}
      // Finally dispose the root transform
      try { hudRoot && hudRoot.dispose(); } catch {}
    }
  };
}
