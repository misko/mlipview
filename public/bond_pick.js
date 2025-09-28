// public/bond_pick.js
// Click-to-select bond + big buttons to rotate via torsion controller.
// Deselect on background click and reattach camera so orbit/pan/zoom work again.

import { makePicker } from "./selection.js";
import { registerBondClear, selectBond, clearSelection as clearGlobalSelection, isAtomDragging } from "./selection_state.js";
import { createEmptySelection, applyBondClick, bondOrientationColor, orientationToSide } from "./selection-model.js";

export function enableBondPicking(
  scene,
  { molecule, torsion, rotatableSpec = [], strict = false }
) {
  // Build whitelist if strict mode enabled
  const allow = new Set();
  if (rotatableSpec && rotatableSpec.length) {
    for (const r of rotatableSpec) {
      const a = Math.min(r.i, r.j), b = Math.max(r.i, r.j);
      allow.add(`${a}-${b}`);
    }
  }
  const isAllowed = (i, j) => !strict || allow.has(`${Math.min(i, j)}-${Math.max(i, j)}`);

  // Selected-bond visual marker (cyan, slightly thicker)
  const selMat = new BABYLON.StandardMaterial("bondSelMat", scene);
  selMat.diffuseColor = new BABYLON.Color3(0.1, 0.9, 1.0);
  selMat.emissiveColor = new BABYLON.Color3(0.05, 0.4, 0.5);

  // Highlight cylinder: 10% larger than base bond diameter (base assumed diameter=1)
  const selMesh = BABYLON.MeshBuilder.CreateCylinder("bondSelMesh", {
    height: 1, diameter: 1.1, tessellation: 32
  }, scene);
  selMesh.material = selMat;
  selMesh.isPickable = false;
  selMesh.setEnabled(false);

  function orientBondMarker(i, j) {
    const ai = molecule.atoms[i].pos, aj = molecule.atoms[j].pos;
    const mid = ai.add(aj).scale(0.5);
    const v = aj.subtract(ai);
    const len = v.length();
    if (len < 1e-6) { selMesh.setEnabled(false); return; }

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

  // Controller-friendly UI
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
    <button id="btnMinus" style="font-size:20px; padding:10px 14px;">⟲ −</button>
    <button id="btnPlus"  style="font-size:20px; padding:10px 14px;">⟳ +</button>
    <!-- recompute button removed -->
  `;
  document.body.appendChild(ui);

  const labelEl = ui.querySelector("#bondLabel");
  // side button removed (cycling handled by repeated selection clicks)
  const btnMinus = ui.querySelector("#btnMinus");
  const btnPlus  = ui.querySelector("#btnPlus");
  // recompute button removed

  const selection = createEmptySelection();
  let selected = null; // local convenience mirror (will hold augmented step/label)

  function updateLabel(warn = false) {
    if (!selected) {
      labelEl.textContent = "Select a bond…";
      // no side button now
      return;
    }
    const { i, j, orientation, step, label } = selected;
    const orientLabel = orientation === 0 ? '(i,j)' : '(j,i)';
    const base = `${label || `bond ${i}-${j}`} ${orientLabel} (step ${step}°)`;
    labelEl.textContent = warn ? `${base} — not in ROY.BONDS` : base;
  }

  function clearSelection() {
    console.log('[bond_pick] clearSelection invoked');
    selected = null;
    selection.kind = null; selection.data = null;
    selMesh.setEnabled(false);
    updateLabel(false);
    // ensure camera reattaches so orbit/pan work immediately
    const cam = scene.activeCamera;
    const canvas = scene.getEngine().getRenderingCanvas();
    if (cam && canvas) cam.attachControl(canvas, true);
  }

  // Register how to clear bond selection so atom selection can invoke it
  registerBondClear(() => {
    clearSelection();
  });

  function applyRotation(sign) {
  if (!selected) return;
  const side = orientationToSide(selected.orientation);
  torsion.rotateAroundBond({ i: selected.i, j: selected.j, orientation: selected.orientation, side, angleDeg: sign * Math.abs(selected.step || 5) });
    try {
    // Snap to float64 authoritative positions and force energy recompute (state simplification)
      if (typeof molecule?.markChanged === 'function') molecule.markChanged();
      if (typeof window.appState?.debugPrint === 'function') window.appState.debugPrint('[desktop rotate]');
    } catch {}
  orientBondMarker(selected.i, selected.j);
  try { if (typeof window.recordEnergyStep === 'function') window.recordEnergyStep(); } catch {}
  }

  btnMinus.onclick = () => applyRotation(-1);
  btnPlus.onclick  = () => applyRotation(+1);
  // side button removed
  // recompute control removed

  // Shortcuts (desktop)
  window.addEventListener("keydown", (e) => {
    if (e.key === "[") applyRotation(-1);
    else if (e.key === "]") applyRotation(+1);
    // side toggle shortcut removed
    else if (e.key === "Escape") clearSelection();
  });

  // Helpers
  const groupByKey = molecule.bondGroups;
  const picker = makePicker(scene, molecule);
  function pickBondAtPointer() {
    const res = picker.pickBondAtPointer();
    if (res) return { key: res.key, mesh: res.mesh, index: res.index };
    return null;
  }
  function isAtomMesh(mesh) { return picker.isAtomMesh(mesh); }

  // Pointer logic:
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERDOWN) return;

  // Do not allow bond selection while atom is being dragged
  if (isAtomDragging && isAtomDragging()) { console.log('[bond_pick] pointerdown ignored: atom dragging'); return; }

  const bondHit = pickBondAtPointer();
    if (bondHit) {
      console.log('[bond_pick] pointerdown bondHit', bondHit.key, 'idx', bondHit.index);
      const group = groupByKey.get(bondHit.key);
      const pair = group?.indices[bondHit.index];
      if (!pair) return;

      const { i, j } = pair;
      // If already selected, advance orientation cycle or clear
      // Use shared selection model
      const prevKind = selection.kind;
      const prevOrientation = selection.data?.orientation;
      const result = applyBondClick(selection, { i, j, key: bondHit.key, index: bondHit.index });
      if (result === 'cleared') {
        console.log('[bond_pick] applyBondClick returned cleared');
        clearSelection();
        return;
      }
      if (selection.kind === 'bond') {
        // augment with step/label spec
        let spec = rotatableSpec.find(r => (r.i === selection.data.i && r.j === selection.data.j) || (r.i === selection.data.j && r.j === selection.data.i));
        if (!spec) spec = { i: selection.data.i, j: selection.data.j, step: 5, label: `bond ${selection.data.i}-${selection.data.j}` };
        selected = { i: selection.data.i, j: selection.data.j, orientation: selection.data.orientation, step: spec.step || 5, label: spec.label };
        const allowed = isAllowed(selected.i, selected.j);
        // set marker color for current orientation
        const col = bondOrientationColor(selected.orientation);
        selMat.diffuseColor = new BABYLON.Color3(col.diffuse.r, col.diffuse.g, col.diffuse.b);
        selMat.emissiveColor = new BABYLON.Color3(col.emissive.r, col.emissive.g, col.emissive.b);
        updateLabel(!allowed);
        orientBondMarker(selected.i, selected.j);
        if (prevKind !== 'bond') { console.log('[bond_pick] selecting bond (prevKind=', prevKind, ')'); selectBond(); }
        return;
      }
    }

    // Not a bond: if it's not an atom either, deselect & reattach camera
    const anyPick = scene.pick(scene.pointerX, scene.pointerY);
    if (!anyPick?.hit || !isAtomMesh(anyPick.pickedMesh)) {
      console.log('[bond_pick] background click clearing selection');
      clearGlobalSelection();
    }
  });
}
