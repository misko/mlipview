// public/interaction.js
// Click/drag atoms with a ground-aligned drag plane. Updates thin instances and bonds.

import { makePicker } from "./selection.js";
import { registerAtomClear, selectAtom, clearSelection as clearGlobalSelection, setAtomDragging } from "./selection_state.js";

export function enableAtomDragging(scene, { atoms, refreshBonds, molecule }) {
  const canvas = scene.getEngine().getRenderingCanvas();
  const camera = scene.activeCamera;

  // --- Debug drag plane (only visible while dragging), never pickable
  const dragPlaneMesh = BABYLON.MeshBuilder.CreateGround("dragPlaneDebug", {
    width: 20, height: 20
  }, scene);
  const dragMat = new BABYLON.StandardMaterial("dragMat", scene);
  dragMat.alpha = 0.12;
  dragMat.diffuseColor = new BABYLON.Color3(0.2, 0.7, 0.9);
  dragPlaneMesh.material = dragMat;
  dragPlaneMesh.setEnabled(false);
  dragPlaneMesh.isPickable = false; // important: don't intercept picks

  // --- Little pick marker (cyan disc)
  const pickMarker = BABYLON.MeshBuilder.CreateDisc("pickMarker", { radius: 0.2, tessellation: 40 }, scene);
  const pmMat = new BABYLON.StandardMaterial("pmMat", scene);
  pmMat.emissiveColor = new BABYLON.Color3(0.1, 0.9, 1.0);
  pmMat.diffuseColor = new BABYLON.Color3(0, 0, 0);
  pickMarker.material = pmMat;
  pickMarker.billboardMode = BABYLON.Mesh.BILLBOARDMODE_Y; // lie flat on plane
  pickMarker.setEnabled(false);
  pickMarker.isPickable = false;

  // --- Atom highlight sphere (cyan glow) shown while dragging
  const atomHighlight = BABYLON.MeshBuilder.CreateSphere("atomPickHighlight", { diameter: 1, segments: 24 }, scene);
  const ahMat = new BABYLON.StandardMaterial("atomPickMat", scene);
  ahMat.diffuseColor = new BABYLON.Color3(0, 0, 0);
  ahMat.emissiveColor = new BABYLON.Color3(0.1, 0.9, 1.0);
  ahMat.alpha = 0.9;
  atomHighlight.material = ahMat;
  atomHighlight.isPickable = false;
  atomHighlight.setEnabled(false);

  // --- Atom hover highlight (subtle)
  const atomHover = BABYLON.MeshBuilder.CreateSphere("atomHoverHighlight", { diameter: 1, segments: 24 }, scene);
  const ahvMat = new BABYLON.StandardMaterial("atomHoverMat", scene);
  ahvMat.diffuseColor = new BABYLON.Color3(0, 0, 0);
  ahvMat.emissiveColor = new BABYLON.Color3(0.06, 0.55, 0.65);
  ahvMat.alpha = 0.6;
  atomHover.material = ahvMat;
  atomHover.isPickable = false;
  atomHover.setEnabled(false);

  // --- Atom persistent selection (cyan ball that stays after click)
  const atomSelect = BABYLON.MeshBuilder.CreateSphere("atomSelect", { diameter: 1, segments: 24 }, scene);
  const asMat = new BABYLON.StandardMaterial("atomSelectMat", scene);
  asMat.diffuseColor = new BABYLON.Color3(0, 0, 0);
  asMat.emissiveColor = new BABYLON.Color3(0.1, 0.9, 1.0);
  asMat.alpha = 0.9;
  atomSelect.material = asMat;
  atomSelect.isPickable = false;
  atomSelect.setEnabled(false);

  // --- Simple tooltip for atom hover
  const atomTip = document.createElement('div');
  atomTip.style.position = 'absolute';
  atomTip.style.pointerEvents = 'none';
  atomTip.style.padding = '2px 6px';
  atomTip.style.borderRadius = '6px';
  atomTip.style.background = 'rgba(15,18,24,0.85)';
  atomTip.style.color = '#cfe9ff';
  atomTip.style.font = '12px system-ui,Segoe UI,Roboto,sans-serif';
  atomTip.style.border = '1px solid rgba(255,255,255,0.08)';
  atomTip.style.transform = 'translate(8px, 8px)';
  atomTip.style.display = 'none';
  atomTip.id = 'atomTooltip';
  document.body.appendChild(atomTip);

  // Ray helpers
  function rayFromPointer() {
    return scene.createPickingRay(scene.pointerX, scene.pointerY, BABYLON.Matrix.Identity(), camera, false);
  }
  function intersectRayPlane(ray, plane) {
    const t = BABYLON.Vector3.Dot(plane.normal, plane.point.subtract(ray.origin)) /
              BABYLON.Vector3.Dot(plane.normal, ray.direction);
    if (!isFinite(t) || t < 0) return null;
    return ray.origin.add(ray.direction.scale(t));
  }

  // Atom picking (shared logic)
  const picker = makePicker(scene, molecule);
  // Helpers: convert atom local position (thin instance center) to world, and get uniform world scale
  function atomLocalToWorld(mesh, localPos) {
    try {
      const wm = mesh.getWorldMatrix();
      const out = new BABYLON.Vector3();
      BABYLON.Vector3.TransformCoordinatesToRef(localPos, wm, out);
      return out;
    } catch {
      return localPos.clone();
    }
  }
  function uniformWorldScale(mesh) {
    try {
      const wm = mesh.getWorldMatrix();
      const s = new BABYLON.Vector3();
      const q = new BABYLON.Quaternion();
      const t = new BABYLON.Vector3();
      wm.decompose(s, q, t);
      return (s.x + s.y + s.z) / 3;
    } catch {
      return 1;
    }
  }
  function pickAtomAtPointer() {
    const res = picker.pickAtomAtPointer();
    if (res) {
      return { atom: res.atom, mesh: res.mesh, index: res.perIndex, globalIndex: res.idx };
    }
    return null;
  }

  // Drag session state
  let dragging = null; // { atom, mesh, index, dragPlane, grabOffset }
  let selected = null; // { atom, mesh, index, globalIndex }
  let pointerDown = null; // { x, y, pick }
  function setCursor(s) { canvas.style.cursor = s; }

  // Register how to clear atom selection so bond UI can invoke it
  registerAtomClear(() => {
    selected = null;
    atomSelect.setEnabled(false);
  });

  function makeDragPlaneThrough(p, normal = new BABYLON.Vector3(0, 1, 0)) {
    return { point: p.clone(), normal: normal.clone() };
  }

  // ----- DOWN: start potential drag -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERDOWN) return;

    const res = pickAtomAtPointer();
    console.log("[pick] DOWN pick:", res);
    if (!res) {
      // Not an atom: do NOT clear selection here; allow bond selection to persist
      return;
    }

  // Record pointer for click/drag threshold
    pointerDown = { x: scene.pointerX, y: scene.pointerY, pick: res };

  // Announce atom selection immediately to clear any bond selection before drag
  selectAtom();

  // Detach camera because we begin a drag session immediately; we will reattach on UP
    if (camera && canvas) camera.detachControl(canvas);

  const plane = makeDragPlaneThrough(res.atom.pos); // horizontal plane through atom
    const ray = rayFromPointer();
    const hit = intersectRayPlane(ray, plane);
    const grabOffset = hit ? hit.subtract(res.atom.pos) : BABYLON.Vector3.Zero();

  dragging = { atom: res.atom, mesh: res.mesh, index: res.index, dragPlane: plane, grabOffset };
  setAtomDragging(true);

    // Show plane + marker
    dragPlaneMesh.position.copyFrom(res.atom.pos);
    dragPlaneMesh.setEnabled(true);
    pickMarker.setEnabled(true);
  // Hide hover highlight/tooltip when we engage drag
  atomHover.setEnabled(false);
  atomTip.style.display = 'none';
  // Show atom highlight (scaled a bit larger than the atom visual size) in world space
  const uw0 = uniformWorldScale(res.mesh);
  const hs = (res.atom.scale || 1) * 1.25 * uw0;
  atomHighlight.scaling.copyFromFloats(hs, hs, hs);
  atomHighlight.position.copyFrom(atomLocalToWorld(res.mesh, res.atom.pos));
    atomHighlight.setEnabled(true);
    setCursor("grabbing");
  });

  // ----- MOVE: update drag -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERMOVE) return;
    // Hover when not dragging
    if (!dragging) {
      const res = picker.pickAtomAtPointer();
      const show = !!res;
      if (show) {
        const uw = uniformWorldScale(res.mesh);
        const hs = (res.atom.scale || 1) * 1.1 * uw;
        atomHover.scaling.copyFromFloats(hs, hs, hs);
        atomHover.position.copyFrom(atomLocalToWorld(res.mesh, res.atom.pos));
        atomHover.setEnabled(true);
        // Tooltip near pointer
        atomTip.textContent = `Atom ${res.type} #${res.idx}`;
        atomTip.style.left = `${scene.pointerX}px`;
        atomTip.style.top = `${scene.pointerY}px`;
        atomTip.style.display = 'block';
      } else {
        atomHover.setEnabled(false);
        atomTip.style.display = 'none';
      }
      return;
    }

    const ray = rayFromPointer();
    const hit = intersectRayPlane(ray, dragging.dragPlane);
    if (!hit) return;

    const targetCenter = hit.subtract(dragging.grabOffset);
    pickMarker.position.copyFrom(targetCenter);
  // Move atom highlight with the dragged atom (world space)
  const uw = uniformWorldScale(dragging.mesh);
  const hs2 = (dragging.atom.scale || 1) * 1.25 * uw;
  atomHighlight.scaling.copyFromFloats(hs2, hs2, hs2);
  atomHighlight.position.copyFrom(atomLocalToWorld(dragging.mesh, targetCenter));

    // Update atom matrix at its per-element index
    const s = dragging.atom.scale ?? 1;
    const mat = BABYLON.Matrix.Compose(
      new BABYLON.Vector3(s, s, s),
      BABYLON.Quaternion.Identity(),
      targetCenter
    );
    molecule.updateAtomMatrixByElement(dragging.atom.type, dragging.index, mat);
    dragging.atom.pos.copyFrom(targetCenter);

    // Fast bond follow (don't mark as changed during drag, only at the end)
    if (typeof refreshBonds === "function") refreshBonds();
  });

  // ----- UP: end drag -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERUP) return;

    if (dragging) {
      // Detect click (no significant movement) vs real drag
      const dx = (pointerDown ? (scene.pointerX - pointerDown.x) : 0);
      const dy = (pointerDown ? (scene.pointerY - pointerDown.y) : 0);
      const clickLike = (dx * dx + dy * dy) <= 9; // <= 3px movement

      const ray = rayFromPointer();
      const hit = intersectRayPlane(ray, dragging.dragPlane);
      if (hit) {
        const finalCenter = hit.subtract(dragging.grabOffset);
        dragging.atom.pos.copyFrom(finalCenter);

        const s = dragging.atom.scale ?? 1;
        const m = BABYLON.Matrix.Compose(
          new BABYLON.Vector3(s, s, s),
          BABYLON.Quaternion.Identity(),
          finalCenter
        );
        molecule.updateAtomMatrixByElement(dragging.atom.type, dragging.index, m);
        
        // Mark molecule as changed for physics cache invalidation
        if (molecule.markChanged) {
          molecule.markChanged();
        }
        
        if (typeof refreshBonds === "function") refreshBonds();

        console.log("[pick] DROP @", finalCenter.toString(), "mesh:", dragging.mesh.name, "index:", dragging.index);
      } else {
        console.log("[pick] DROP: no plane hit (keeping last position)");
      }

      // If it was effectively a click, persist selection on this atom
      if (pointerDown && clickLike) {
        selected = { atom: dragging.atom, mesh: dragging.mesh, index: dragging.index, globalIndex: pointerDown.pick.globalIndex };
        const uw = uniformWorldScale(selected.mesh);
        const hs = (selected.atom.scale || 1) * 1.2 * uw;
        atomSelect.scaling.copyFromFloats(hs, hs, hs);
        atomSelect.position.copyFrom(atomLocalToWorld(selected.mesh, selected.atom.pos));
        atomSelect.setEnabled(true);
        // Announce atom selection to clear any bond selection
        selectAtom();
      } else {
        // After a real drag, keep the atom selected as well
        selected = { atom: dragging.atom, mesh: dragging.mesh, index: dragging.index, globalIndex: pointerDown?.pick?.globalIndex };
        const uw = uniformWorldScale(selected.mesh);
        const hs = (selected.atom.scale || 1) * 1.2 * uw;
        atomSelect.scaling.copyFromFloats(hs, hs, hs);
        atomSelect.position.copyFrom(atomLocalToWorld(selected.mesh, selected.atom.pos));
        atomSelect.setEnabled(true);
        // Announce atom selection to clear any bond selection
        selectAtom();
      }

  dragging = null;
  setAtomDragging(false);
      pickMarker.setEnabled(false);
      dragPlaneMesh.setEnabled(false);
      atomHighlight.setEnabled(false);
      atomHover.setEnabled(false);
      atomTip.style.display = 'none';
      setCursor("default");
      if (camera && canvas) camera.attachControl(canvas, true);
      console.log("[pick] Camera reattached");
      pointerDown = null;
    } else {
      // Safety: ensure camera is attached even if we weren't dragging
      if (camera && canvas) camera.attachControl(canvas, true);
      // If this was a true background click (neither atom nor bond), clear selection
      const anyPick = scene.pick(scene.pointerX, scene.pointerY);
      const mesh = anyPick?.hit ? anyPick.pickedMesh : null;
      const isAtom = mesh && picker.isAtomMesh(mesh);
      const isBond = mesh && picker.isBondMesh(mesh);
      if (!mesh || (!isAtom && !isBond)) {
        clearGlobalSelection();
      }
    }
  });

  // Keep persistent selection synced with atom/world transforms
  scene.onBeforeRenderObservable.add(() => {
    if (!selected) return;
    try {
      const uw = uniformWorldScale(selected.mesh);
      const hs = (selected.atom.scale || 1) * 1.2 * uw;
      atomSelect.scaling.copyFromFloats(hs, hs, hs);
      atomSelect.position.copyFrom(atomLocalToWorld(selected.mesh, selected.atom.pos));
      if (!atomSelect.isEnabled()) atomSelect.setEnabled(true);
    } catch {
      // noop
    }
  });
}
// (no-op: helper duplicates removed)
