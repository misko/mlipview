// public/bond_pick.js
// Click-to-select bond + big buttons to rotate via torsion controller.
// Now allows any bond by default (spec is only enforced if strict=true).

export function enableBondPicking(
  scene,
  { molecule, torsion, rotatableSpec = [], strict = false } // <-- strict OFF by default
) {
  const canvas = scene.getEngine().getRenderingCanvas();

  // Build a quick "allowed" set if a spec is present
  const allow = new Set();
  if (rotatableSpec && rotatableSpec.length) {
    for (const r of rotatableSpec) {
      const a = Math.min(r.i, r.j), b = Math.max(r.i, r.j);
      allow.add(`${a}-${b}`);
    }
  }

  function isAllowed(i, j) {
    if (!strict) return true; // default: allow any bond
    const a = Math.min(i, j), b = Math.max(i, j);
    return allow.has(`${a}-${b}`);
  }

  // Visual marker for selected bond (slightly thicker cyan cylinder)
  const selMat = new BABYLON.StandardMaterial("bondSelMat", scene);
  selMat.diffuseColor = new BABYLON.Color3(0.1, 0.9, 1.0);
  selMat.emissiveColor = new BABYLON.Color3(0.05, 0.4, 0.5);

  const selMesh = BABYLON.MeshBuilder.CreateCylinder("bondSelMesh", {
    height: 1, diameter: 1.25, tessellation: 32
  }, scene);
  selMesh.material = selMat;
  selMesh.isPickable = false;
  selMesh.setEnabled(false);

  function orientBondMarker(i, j) {
    const ai = molecule.atoms[i].pos;
    const aj = molecule.atoms[j].pos;
    const mid = ai.add(aj).scale(0.5);
    const v = aj.subtract(ai);
    const len = v.length();
    if (len < 1e-6) return selMesh.setEnabled(false);

    const up = BABYLON.Vector3.Up();
    const d = v.normalizeToNew();
    const dot = BABYLON.Vector3.Dot(up, d);
    let rot;
    if (dot > 0.9999) rot = BABYLON.Quaternion.Identity();
    else if (dot < -0.9999) rot = BABYLON.Quaternion.RotationAxis(BABYLON.Vector3.Right(), Math.PI);
    else {
      const axis = BABYLON.Vector3.Cross(up, d).normalize();
      rot = BABYLON.Quaternion.RotationAxis(axis, Math.acos(dot));
    }
    selMesh.scaling.setAll(1);
    selMesh.scaling.y = len;
    selMesh.position.copyFrom(mid);
    selMesh.rotationQuaternion = rot;
    selMesh.setEnabled(true);
  }

  // UI: big VR-friendly buttons
  const ui = document.createElement("div");
  ui.style.position = "absolute";
  ui.style.bottom = "12px";
  ui.style.left = "50%";
  ui.style.transform = "translateX(-50%)";
  ui.style.background = "rgba(15,18,24,0.75)";
  ui.style.border = "1px solid rgba(255,255,255,0.08)";
  ui.style.borderRadius = "10px";
  ui.style.padding = "10px 12px";
  ui.style.display = "flex";
  ui.style.gap = "8px";
  ui.style.alignItems = "center";
  ui.style.color = "#d7e6ff";
  ui.style.font = "14px system-ui,-apple-system,Segoe UI,Roboto,sans-serif";
  ui.style.userSelect = "none";
  ui.style.pointerEvents = "auto";
  ui.style.backdropFilter = "blur(4px)";
  ui.innerHTML = `
    <span id="bondLabel" style="min-width: 260px; display:inline-block;">Select a bond…</span>
    <button id="btnSide" style="font-size:16px; padding:10px 12px;">Side: j</button>
    <button id="btnMinus" style="font-size:20px; padding:10px 14px;">⟲ −</button>
    <button id="btnPlus"  style="font-size:20px; padding:10px 14px;">⟳ +</button>
    <button id="btnRecompute" title="Recompute connectivity after next rotation">recompute</button>
  `;
  document.body.appendChild(ui);

  const labelEl = ui.querySelector("#bondLabel");
  const btnSide = ui.querySelector("#btnSide");
  const btnMinus = ui.querySelector("#btnMinus");
  const btnPlus  = ui.querySelector("#btnPlus");
  const btnRec   = ui.querySelector("#btnRecompute");

  let selected = null; // { i, j, side, step, label }
  let doRecompute = false;

  function updateLabel(warnNotInSpec = false) {
    if (!selected) {
      labelEl.textContent = "Select a bond…";
      btnSide.textContent = "Side: j";
      return;
    }
    const { i, j, side, step, label } = selected;
    const base = `${label || `bond ${i}-${j}`} (step ${step}°)`;
    labelEl.textContent = warnNotInSpec ? `${base} — not in ROY.BONDS` : base;
    btnSide.textContent = `Side: ${side}`;
  }

  function applyRotation(sign) {
    if (!selected) return;
    torsion.rotateAroundBond({
      i: selected.i,
      j: selected.j,
      side: selected.side,
      angleDeg: sign * Math.abs(selected.step || 5),
      recompute: doRecompute
    });
    doRecompute = false;
    orientBondMarker(selected.i, selected.j);
  }

  btnMinus.onclick = () => applyRotation(-1);
  btnPlus.onclick  = () => applyRotation(+1);
  btnSide.onclick  = () => {
    if (!selected) return;
    selected.side = (selected.side === "j") ? "i" : "j";
    updateLabel();
  };
  btnRec.onclick = () => { doRecompute = true; };

  // Keyboard shortcuts: [ / ] / S
  window.addEventListener("keydown", (e) => {
    if (e.key === "[") applyRotation(-1);
    else if (e.key === "]") applyRotation(+1);
    else if (e.key.toLowerCase() === "s") btnSide.click();
  });

  function pickBondAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY, (m) => m.name?.startsWith("bond_"));
    if (pick?.hit && pick.pickedMesh && pick.thinInstanceIndex != null && pick.thinInstanceIndex >= 0) {
      const mesh = pick.pickedMesh;
      const key  = mesh.name.replace(/^bond_/, "");
      const idx  = pick.thinInstanceIndex;
      return { key, mesh, index: idx };
    }
    return null;
  }

  const groupByKey = molecule.bondGroups;

  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERDOWN) return;
    const res = pickBondAtPointer();
    if (!res) return;

    const group = groupByKey.get(res.key);
    if (!group) return;
    const pair = group.indices[res.index];
    if (!pair) return;

    const { i, j } = pair;

    // Log indices to help authoring ROY.BONDS
    console.log(`[bond-pick] Selected bond i=${i} j=${j} (pairKey=${res.key})`);

    const allowed = isAllowed(i, j);
    if (!allowed) {
      // Show marker but warn; still select (since strict=true would block rotation anyway)
      orientBondMarker(i, j);
      // Use a default step and side so user can flip side if needed
      selected = { i, j, side: "j", step: 5, label: `bond ${i}-${j}` };
      updateLabel(true);
      return;
    }

    // If spec has a matching entry, use it; else sensible defaults
    let spec = rotatableSpec.find(r => (r.i === i && r.j === j) || (r.i === j && r.j === i));
    if (!spec) spec = { i, j, side: "j", step: 5, label: `bond ${i}-${j}` };

    selected = { i: spec.i, j: spec.j, side: spec.side || "j", step: spec.step || 5, label: spec.label };
    updateLabel(false);
    orientBondMarker(i, j);
  });
}
