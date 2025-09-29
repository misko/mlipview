// mlipviewer2/public/interaction/adapters/desktopInputAdapter.js
// Translates mouse events into dragCore + selection actions.

export function createDesktopInputAdapter({ scene, canvas, dragCore, pickingFacade, selectionService }) {
  let pointerDownInfo = null; // { x,y, atomIndex }
  const CLICK_EPS2 = 9; // px^2 threshold

  function onPointerDown(evt) {
    const atom = pickingFacade.pickAtomAtPointer();
    const bond = !atom ? pickingFacade.pickBondAtPointer() : null;
    if (atom) {
      pointerDownInfo = { x: scene.pointerX, y: scene.pointerY, atomIndex: atom.index };
      dragCore.begin(atom.index);
    } else if (bond) {
      selectionService.selectBond({ i: bond.i, j: bond.j, orientation: 0 });
    } else {
      selectionService.clearSelection();
    }
  }

  function onPointerMove(evt) {
    if (dragCore.isDragging()) {
      // Project pointer to plane through starting atom's current y (simple horizontal plane); could improve later.
      const ray = scene.createPickingRay(scene.pointerX, scene.pointerY, BABYLON.Matrix.Identity(), scene.activeCamera);
      // Horizontal plane y = current atom y
      const sel = selectionService.getSelection();
      if (!sel || sel.kind !== 'atom') return;
      const idx = sel.data.index;
      const baseY = dragCore.isDragging() ? 0 + (/* future: store original y */ 0) : 0;
      const planeNormal = new BABYLON.Vector3(0,1,0);
      const denom = BABYLON.Vector3.Dot(planeNormal, ray.direction);
      if (Math.abs(denom) < 1e-6) return;
      const t = -(BABYLON.Vector3.Dot(planeNormal, ray.origin) - baseY) / denom;
      if (t < 0) return;
      const hit = ray.origin.add(ray.direction.scale(t));
      dragCore.move({ x: hit.x, y: hit.y, z: hit.z });
    }
  }

  function onPointerUp(evt) {
    if (dragCore.isDragging()) {
      const dx = scene.pointerX - pointerDownInfo.x;
      const dy = scene.pointerY - pointerDownInfo.y;
      const clickLike = (dx*dx + dy*dy) <= CLICK_EPS2;
      dragCore.end(true);
      if (clickLike) {
        // Already selected as atom; nothing extra
      }
    }
    pointerDownInfo = null;
  }

  scene.onPointerObservable.add((pi) => {
    if (pi.type === BABYLON.PointerEventTypes.POINTERDOWN) onPointerDown(pi.event);
    else if (pi.type === BABYLON.PointerEventTypes.POINTERMOVE) onPointerMove(pi.event);
    else if (pi.type === BABYLON.PointerEventTypes.POINTERUP) onPointerUp(pi.event);
  });

  return { dispose() {/* future: remove observers */} };
}
