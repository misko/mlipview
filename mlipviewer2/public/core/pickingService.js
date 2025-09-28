// Picking & Drag Service
// Responsibilities:
//  1. Pointer picking (atoms/bonds) and selection updates.
//  2. Initiate atom drag when an atom is selected (desktop & VR can both supply an intersector pattern).
//  3. Temporarily freeze ArcRotateCamera during active drag by detaching controls and zeroing inertial offsets each frame.
//  4. Reattach camera controls on drag end.
// Notes for VR integration:
//  - VR controller ray logic should call the same selection methods before invoking manipulation.beginDrag.
//  - Intersector currently uses a simple pointer->plane projection; VR layer can override by calling manipulation.setDragPlane then updateDrag with its own ray-plane hits.
//  - This module intentionally avoids any DOM-specific APIs except canvas event listeners so it can be adapted for XR sessions.
// Options param: { manipulation, camera }
export function createPickingService(scene, view, selectionService, { manipulation, camera } = {}) {
  let cameraDetachedForDrag = false;
  let cameraLock = null; // { alpha,beta,radius,target }
  // Helper to fully freeze camera even if Babylon processes some inputs (inertia or event ordering)
  function freezeCameraFrame() {
    if (!cameraDetachedForDrag || !camera) return;
    // Force inertial offsets to zero to stop momentum updates
    if ('inertialAlphaOffset' in camera) camera.inertialAlphaOffset = 0;
    if ('inertialBetaOffset' in camera) camera.inertialBetaOffset = 0;
    if ('inertialRadiusOffset' in camera) camera.inertialRadiusOffset = 0;
    if ('inertialPanningX' in camera) camera.inertialPanningX = 0;
    if ('inertialPanningY' in camera) camera.inertialPanningY = 0;
    // Hard lock: restore original params in case any sneak through
    if (cameraLock) {
      if (camera.alpha !== cameraLock.alpha) camera.alpha = cameraLock.alpha;
      if (camera.beta !== cameraLock.beta) camera.beta = cameraLock.beta;
      if (camera.radius !== cameraLock.radius) camera.radius = cameraLock.radius;
      if (camera.target && cameraLock.target) {
        if (camera.target.x !== cameraLock.target.x) camera.target.x = cameraLock.target.x;
        if (camera.target.y !== cameraLock.target.y) camera.target.y = cameraLock.target.y;
        if (camera.target.z !== cameraLock.target.z) camera.target.z = cameraLock.target.z;
      }
    }
  }
  function pickAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY);
    if (!pick || !pick.hit) return null;
    const atom = view.resolveAtomPick(pick);
    if (atom) return atom;
    const bond = view.resolveBondPick(pick);
    if (bond) return bond;
    return null;
  }
  function handlePointerDown() {
    const res = pickAtPointer();
    if (!res) { selectionService.clear(); return; }
  if (res.kind === 'atom') selectionService.clickAtom(res.index);
    else if (res.kind === 'bond') selectionService.clickBond(res);
    // After selection, if atom selected and manipulation present, initiate drag baseline
    if (manipulation && selectionService.get().kind === 'atom') {
      // Provide an intersector that projects pointer ray onto Z=0 plane for now
      const intersector = (planePoint, planeNormal) => {
        // Build a simple ray from camera via scene.unproject if available; fallback: x,y with fixed z
        // Simplified: treat pointer coords as x,y in plane units (debug only). Real implementation would use actual picking ray.
        return { x: scene.pointerX/50, y: -scene.pointerY/50, z: 0 };
      };
      const started = manipulation.beginDrag(intersector);
      if (started && camera && camera.detachControl && !cameraDetachedForDrag) {
        const canvas = scene.getEngine && scene.getEngine().getRenderingCanvas ? scene.getEngine().getRenderingCanvas() : undefined;
        try { camera.detachControl(canvas); cameraDetachedForDrag = true; } catch (e) { /* ignore */ }
        cameraLock = { alpha: camera.alpha, beta: camera.beta, radius: camera.radius, target: camera.target && { x: camera.target.x, y: camera.target.y, z: camera.target.z } };
      }
    }
  }
  function handlePointerMove() {
    if (!manipulation || !cameraDetachedForDrag) return;
    const intersector = () => ({ x: scene.pointerX/50, y: -scene.pointerY/50, z:0 });
    manipulation.updateDrag(intersector);
  }
  function handlePointerUp() {
    if (!manipulation || !cameraDetachedForDrag) return;
    manipulation.endDrag();
    if (camera && camera.attachControl) {
      const canvas = scene.getEngine && scene.getEngine().getRenderingCanvas ? scene.getEngine().getRenderingCanvas() : undefined;
      try { camera.attachControl(canvas, true); } catch (e) { /* ignore */ }
    }
    cameraDetachedForDrag = false;
    cameraLock = null;
  }
  scene.onPointerObservable.add(pi => {
    if (pi.type === BABYLON.PointerEventTypes.POINTERDOWN) handlePointerDown();
  });
  // Attach raw DOM events if available for move/up (browser context). In test stubs these may be absent.
  if (scene.getEngine && scene.getEngine().getRenderingCanvas) {
    const canvas = scene.getEngine().getRenderingCanvas();
    canvas && canvas.addEventListener('pointermove', handlePointerMove);
    canvas && canvas.addEventListener('pointerup', handlePointerUp);
    canvas && canvas.addEventListener('pointerleave', handlePointerUp);
  }
  // Per-frame observer to enforce freeze while dragging
  scene.onBeforeRenderObservable && scene.onBeforeRenderObservable.add(freezeCameraFrame);
  return { pickAtPointer };
}
