// public/interaction.js
// Enables hover + click/drag movement of atoms.
// Requires: enableAtomDragging(scene, { atoms, refreshBonds, molecule })

export function enableAtomDragging(scene, { atoms, refreshBonds, molecule }) {
  const camera = scene.activeCamera;
  const canvas = scene.getEngine().getRenderingCanvas();

  // ----- Visual markers (no mesh-wide highlight) -----
  // Hover marker: thin yellow torus
  const hoverMarker = BABYLON.MeshBuilder.CreateTorus("hoverMarker", {
    diameter: 0.8, thickness: 0.03, tessellation: 24
  }, scene);
  const hoverMat = new BABYLON.StandardMaterial("hoverMat", scene);
  hoverMat.diffuseColor = new BABYLON.Color3(1, 1, 0);
  hoverMat.emissiveColor = new BABYLON.Color3(0.3, 0.3, 0.0);
  hoverMarker.material = hoverMat;
  hoverMarker.setEnabled(false);

  // Picked marker: cyan torus
  const pickMarker = BABYLON.MeshBuilder.CreateTorus("pickMarker", {
    diameter: 0.9, thickness: 0.05, tessellation: 24
  }, scene);
  const pickMat = new BABYLON.StandardMaterial("pickMat", scene);
  pickMat.diffuseColor = new BABYLON.Color3(0, 1, 1);
  pickMat.emissiveColor = new BABYLON.Color3(0.0, 0.3, 0.3);
  pickMarker.material = pickMat;
  pickMarker.setEnabled(false);

  // Visible drag plane for debugging (aligned to camera on pointer down)
  const dragPlaneMesh = BABYLON.MeshBuilder.CreateGround("dragPlaneDebug", {
    width: 20, height: 20
  }, scene);
  const dragMat = new BABYLON.StandardMaterial("dragMat", scene);
  dragMat.alpha = 0.12;
  dragMat.diffuseColor = new BABYLON.Color3(0.2, 0.7, 0.9);
  dragPlaneMesh.material = dragMat;
  dragPlaneMesh.setEnabled(false);

  // ----- Lookup: (mesh.uniqueId, per-mesh index) -> atom -----
  const atomByMeshIndex = new Map();
  for (const a of atoms) atomByMeshIndex.set(`${a.mesh.uniqueId}:${a.index}`, a);

  // ----- State -----
  let dragging = null; // { mesh, index, atom, grabOffset, dragPlane:{n,d} }
  let hoverKey = null;

  // ----- Helpers -----
  function setCursor(c) { if (canvas) canvas.style.cursor = c; }
  function log(...args) { console.log("[pick]", ...args); }

  function makeDragPlane(normal, point) {
    const n = normal.normalizeToNew();
    const d = -BABYLON.Vector3.Dot(n, point);
    return { n, d };
  }
  function rayFromPointer() {
    return scene.createPickingRay(
      scene.pointerX, scene.pointerY, BABYLON.Matrix.Identity(),
      camera, false, camera.viewport
    );
  }
  function intersectRayPlane(ray, plane) {
    const denom = BABYLON.Vector3.Dot(plane.n, ray.direction);
    if (Math.abs(denom) < 1e-6) return null;
    const t = -(BABYLON.Vector3.Dot(plane.n, ray.origin) + plane.d) / denom;
    if (t < 0) return null;
    return ray.origin.add(ray.direction.scale(t));
  }

  // Prefer thin-instance pick; fallback to master (choose nearest instance)
  function pickAtomAtPointer() {
    // 1) Thin instances
    let pick = scene.pick(
      scene.pointerX, scene.pointerY,
      m => m.thinInstanceEnablePicking === true
    );
    if (pick?.hit && pick.pickedMesh && pick.thinInstanceIndex != null && pick.thinInstanceIndex >= 0) {
      return { pick, mesh: pick.pickedMesh, index: pick.thinInstanceIndex, mode: "thin" };
    }
    // 2) Master fallback
    pick = scene.pick(
      scene.pointerX, scene.pointerY,
      m => m.name?.startsWith("base_")
    );
    if (pick?.hit && pick.pickedMesh) {
      const mesh = pick.pickedMesh;
      const candidates = atoms.filter(a => a.mesh === mesh).map(a => a.index);
      let best = { i: -1, d2: Infinity };
      for (const idx of candidates) {
        const key = `${mesh.uniqueId}:${idx}`;
        const a = atomByMeshIndex.get(key);
        if (!a) continue;
        const d2 = BABYLON.Vector3.DistanceSquared(a.pos, pick.pickedPoint);
        if (d2 < best.d2) best = { i: idx, d2 };
      }
      if (best.i >= 0) return { pick, mesh, index: best.i, mode: "master-fallback" };
    }
    return null;
    }

  // ----- HOVER -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERMOVE || dragging) return;
    const res = pickAtomAtPointer();
    if (res) {
      const key = `${res.mesh.uniqueId}:${res.index}`;
      if (hoverKey !== key) {
        hoverKey = key;
        const atom = atomByMeshIndex.get(key);
        if (atom) {
          hoverMarker.position.copyFrom(atom.pos);
          hoverMarker.setEnabled(true);
          setCursor("pointer");
        }
      }
    } else {
      hoverKey = null;
      hoverMarker.setEnabled(false);
      setCursor("default");
    }
  });

  // ----- DOWN: begin drag -----
  scene.onPointerObservable.add((pi) => {
    if (pi.type !== BABYLON.PointerEventTypes.POINTERDOWN) return;
    const res = pickAtomAtPointer();
    log("DOWN pick:", res);
    if (!res) return;

    camera.detachControl(canvas); // stop camera during drag

    const { mesh, index, pick } = res;
    const key = `${mesh.uniqueId}:${index}`;
    const atom = atomByMeshIndex.get(key);
    if (!atom) { log("No atom for", key); return; }

    // Drag plane through atom center, facing camera
    const camForward = camera.getForwardRay().direction;
    const dragPlane = makeDragPlane(camForward, atom.pos);
    const hitPoint = pick?.pickedPoint ?? atom.pos;
    const grabOffset = hitPoint.subtract(atom.pos);

    dragging = { mesh, index, atom, grabOffset, dragPlane };

    // Visuals
    pickMarker.position.copyFrom(atom.pos);
    pickMarker.setEnabled(true);

    const up = new BABYLON.Vector3(0, 1, 0);
    const right = BABYLON.Vector3.Cross(camForward, up).normalize();
    const camOrthoUp = BABYLON.Vector3.Cross(right, camForward).normalize();
    dragPlaneMesh.position.copyFrom(atom.pos);
    dragPlaneMesh.rotation = BABYLON.Vector3.RotationFromAxis(right, camOrthoUp, camForward);
    dragPlaneMesh.setEnabled(true);

    setCursor("grabbing");
  });

  // ----- MOVE: drag -----
  scene.onPointerObservable.add((pi) => {
    if (!dragging || pi.type !== BABYLON.PointerEventTypes.POINTERMOVE) return;

    const ray = rayFromPointer();
    const hit = intersectRayPlane(ray, dragging.dragPlane);
    if (!hit) return;

    const newCenter = hit.subtract(dragging.grabOffset);
    dragging.atom.pos.copyFrom(newCenter);

    // Preserve per-atom scale and bulk-update the element's buffer for rock-solid visuals
    const s = dragging.atom.scale ?? 1;
    const m = BABYLON.Matrix.Compose(
      new BABYLON.Vector3(s, s, s),
      BABYLON.Quaternion.Identity(),
      newCenter
    );
    molecule.updateAtomMatrixByElement(dragging.atom.type, dragging.index, m);

    // Recompute bonds
    refreshBonds();

    // Move markers
    pickMarker.position.copyFrom(newCenter);
    dragPlaneMesh.position.copyFrom(newCenter);
  });

  // ----- UP: end drag (final commit) -----
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

        refreshBonds();
		molecule.recomputeBonds();
        console.log("[pick] DROP @", finalCenter.toString(), "mesh:", dragging.mesh.name, "index:", dragging.index);
      } else {
        console.log("[pick] DROP: no plane hit (keeping last position)");
      }

      // Clear visuals & restore camera
      dragging = null;
      pickMarker.setEnabled(false);
      dragPlaneMesh.setEnabled(false);
      setCursor("default");
      camera.attachControl(canvas, true);
      console.log("[pick] Camera reattached");
    }
  });
}
