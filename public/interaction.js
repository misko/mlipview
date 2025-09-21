// public/interaction.js
// Click/drag atoms with a ground-aligned drag plane. Updates thin instances and bonds.

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

  // Atom picking
  function pickAtomAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY, (m) => m.name?.startsWith("base_"));
    if (pick?.hit && pick.pickedMesh && pick.thinInstanceIndex != null && pick.thinInstanceIndex >= 0) {
      const mesh = pick.pickedMesh; // the element master
      const instanceIndex = pick.thinInstanceIndex;
      // Find which atom this per-mesh index corresponds to (we stored it in molecule.atoms entries)
      const element = mesh.name.replace(/^base_/, "");
      // locate the atom in atoms[] with matching mesh/type/index
      for (let k = 0; k < atoms.length; k++) {
        const a = atoms[k];
        if (a.type === element && a.mesh === mesh && a.index === instanceIndex) {
          return { atom: a, mesh, index: instanceIndex, globalIndex: k };
        }
      }
    }
    return null;
  }

  // Drag session state
  let dragging = null; // { atom, mesh, index, dragPlane, grabOffset }
  function setCursor(s) { canvas.style.cursor = s; }

  function makeDragPlaneThrough(p, normal = new BABYLON.Vector3(0, 1, 0)) {
    return { point: p.clone(), normal: normal.clone() };
  }

  // ----- DOWN: start potential drag -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERDOWN) return;

    const res = pickAtomAtPointer();
    console.log("[pick] DOWN pick:", res);
    if (!res) {
      // Background or non-atom: do NOT detach camera here
      return;
    }

    // Detach camera only when actually dragging an atom
    if (camera && canvas) camera.detachControl(canvas);

    const plane = makeDragPlaneThrough(res.atom.pos); // horizontal plane through atom
    const ray = rayFromPointer();
    const hit = intersectRayPlane(ray, plane);
    const grabOffset = hit ? hit.subtract(res.atom.pos) : BABYLON.Vector3.Zero();

    dragging = { atom: res.atom, mesh: res.mesh, index: res.index, dragPlane: plane, grabOffset };

    // Show plane + marker
    dragPlaneMesh.position.copyFrom(res.atom.pos);
    dragPlaneMesh.setEnabled(true);
    pickMarker.setEnabled(true);
    setCursor("grabbing");
  });

  // ----- MOVE: update drag -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERMOVE) return;
    if (!dragging) return;

    const ray = rayFromPointer();
    const hit = intersectRayPlane(ray, dragging.dragPlane);
    if (!hit) return;

    const targetCenter = hit.subtract(dragging.grabOffset);
    pickMarker.position.copyFrom(targetCenter);

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

      dragging = null;
      pickMarker.setEnabled(false);
      dragPlaneMesh.setEnabled(false);
      setCursor("default");
      if (camera && canvas) camera.attachControl(canvas, true);
      console.log("[pick] Camera reattached");
    } else {
      // Safety: ensure camera is attached even if we weren't dragging
      if (camera && canvas) camera.attachControl(canvas, true);
    }
  });
}
