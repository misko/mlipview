import { setThinInstanceMatrix } from "./helpers.js";

// Simple drag controller that works with mouse/touch.
// In VR, Babylon's Pointer Selection also triggers pointer events in most cases.
// If needed, we can hook WebXR feature observables later.
export function enableAtomDragging(scene, {
  atoms,              // [{ mesh, type, index, pos: Vector3 }]
  refreshBonds        // function()
}) {
  const camera = scene.activeCamera;

  // Build a pickable list (only atom masters; picking returns thinInstanceIndex)
  const pickTargets = new Set(atoms.map(a => a.mesh));

  // State
  let dragging = null;    // { mesh, index, atom, grabOffset, dragPlane }
  const tmp = {
    pickPoint: new BABYLON.Vector3(),
    plane: new BABYLON.Plane(0,1,0,0),
    ray: new BABYLON.Ray(),
    v: new BABYLON.Vector3(),
  };

  // Create a reusable invisible drag plane mesh (we'll compute intersections analytically instead)
  function makeDragPlane(normal, point) {
    // Plane: n·x + d = 0  with d = -n·p0
    const n = normal.normalizeToNew();
    const d = -BABYLON.Vector3.Dot(n, point);
    return { n, d };
  }

  function rayFromPointer(ev) {
    // Build a picking ray from the active camera through the pointer
    return scene.createPickingRay(scene.pointerX, scene.pointerY, BABYLON.Matrix.Identity(), camera, false, camera.viewport);
  }

  function intersectRayPlane(ray, plane) {
    // Solve for t: ray.origin + t*ray.direction hits plane n·x + d = 0
    const denom = BABYLON.Vector3.Dot(plane.n, ray.direction);
    if (Math.abs(denom) < 1e-6) return null;
    const t = -(BABYLON.Vector3.Dot(plane.n, ray.origin) + plane.d) / denom;
    if (t < 0) return null;
    return ray.origin.add(ray.direction.scale(t));
  }

  // Pointer down: pick an atom thin instance
  scene.onPointerObservable.add((pointerInfo) => {
    if (pointerInfo.type !== BABYLON.PointerEventTypes.POINTERDOWN) return;

    const pick = scene.pick(scene.pointerX, scene.pointerY, (mesh) => pickTargets.has(mesh));
    if (!pick?.hit || pick.thinInstanceIndex == null || pick.thinInstanceIndex < 0) return;

    const mesh = pick.pickedMesh;
    const idx = pick.thinInstanceIndex;

    // Find the logical atom (we know which mesh & index)
    const atom = atoms.find(a => a.mesh === mesh && a.index === idx);
    if (!atom) return;

    // Create a drag plane facing the camera, going through the hit point
    const camForward = camera.getForwardRay().direction;
    const dragPlane = makeDragPlane(camForward, pick.pickedPoint);

    // Store grab offset (local offset from atom center to hit point)
    const grabOffset = pick.pickedPoint.subtract(atom.pos);

    dragging = { mesh, index: idx, atom, grabOffset, dragPlane };
  });

  // Pointer move: update atom position
  scene.onPointerObservable.add((pointerInfo) => {
    if (!dragging || pointerInfo.type !== BABYLON.PointerEventTypes.POINTERMOVE) return;

    const ray = rayFromPointer(pointerInfo);
    const hit = intersectRayPlane(ray, dragging.dragPlane);
    if (!hit) return;

    // Maintain the grab offset so it feels natural
    const newCenter = hit.subtract(dragging.grabOffset);

    // Update JS position (this drives bonds too)
    dragging.atom.pos.copyFrom(newCenter);

    // Compose instance matrix (unit scale, identity rot)
    const m = BABYLON.Matrix.Compose(BABYLON.Vector3.One(), BABYLON.Quaternion.Identity(), newCenter);
    setThinInstanceMatrix(dragging.mesh, dragging.index, m);

    // Update connected bonds
    refreshBonds();
  });

  // Pointer up: release
  scene.onPointerObservable.add((pointerInfo) => {
    if (pointerInfo.type !== BABYLON.PointerEventTypes.POINTERUP) return;
    dragging = null;
  });

  // OPTIONAL: enable Babylon gizmos later if you want axis-constrained drags

  // WEBXR NOTE:
  // With the default XR experience + PointerSelection, the controller rays drive scene picks,
  // so this same handler works in VR on Quest (trigger acts like click).
  // If you find a browser that doesn't route those to pointer observables,
  // we can hook into XR's PointerSelection feature observables as a follow-up.
}
