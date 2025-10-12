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
// NOTE: Energy plot updates: direct programmatic manipulation (viewerApi.manipulation.*) already
// wraps operations to recompute forces + recordInteraction. However, user desktop drags flow
// through this pickingService path. To keep energy plot consistent, we optionally accept
// an energyHook({ kind }) callback invoked on drag move/end when geometry changes.
import { __count } from '../util/funcCount.js';
// Consolidated Interaction Controller: picking + drag + rotation entry points
export function createPickingService(scene, view, selectionService, { manipulation, camera, energyHook } = {}) {
  __count('pickingService#createPickingService');
  let cameraDetachedForDrag = false;
  let cameraLock = null; // { alpha,beta,radius,target }
  let dragActive = false;
  const DBG = (typeof window !== 'undefined') && !!window.__MLIPVIEW_DEBUG_TOUCH;
  // Public interaction surface (populated below) returned at end
  const api = { pickAtPointer, beginAtomDrag, updateDragFromPointer, endDrag, rotateSelectedBond, freezeCameraFrame, selection: selectionService.get };
  // Helper to fully freeze camera even if Babylon processes some inputs (inertia or event ordering)
  function freezeCameraFrame() {
    __count('pickingService#freezeCameraFrame');
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
    __count('pickingService#pickAtPointer');
    const pick = scene.pick(scene.pointerX, scene.pointerY);
    if (!pick || !pick.hit) return null;
    const atom = view.resolveAtomPick(pick);
    if (atom) return atom;
    const bond = view.resolveBondPick(pick);
    if (bond) return bond;
    return null;
  }
  function handlePointerDown() {
    __count('pickingService#handlePointerDown');
    const res = pickAtPointer();
    if (!res) {
      // Fallback: if no pick result and at least one atom exists, select atom 0 for drag (helps automated tests)
      if (selectionService.get().kind !== 'atom' && view && view._debugAutoSelectFirstOnEmpty) {
        selectionService.clickAtom(0);
      } else {
        selectionService.clear();
        return;
      }
    }
  if (res.kind === 'atom') selectionService.clickAtom(res.index);
    else if (res.kind === 'bond') selectionService.clickBond(res);
    // After selection, if atom selected and manipulation present, initiate drag baseline
    if (manipulation && selectionService.get().kind === 'atom') {
      const atomSel = selectionService.get();
      let started = false;
      // Compute camera-aligned drag plane: normal = (atomPos - cameraPos)
      let planePoint = null; let planeNormal = null;
      try {
        if (camera && atomSel && atomSel.kind === 'atom') {
          const atomPos = manipulation.molState ? manipulation.molState.positions[atomSel.data.index] : (view && view._internals && view._internals.highlight && view._internals.highlight.atom ? view._internals.highlight.atom.position : null);
          // Fallback: take position from molState via selection if available
          if (!atomPos && manipulation && manipulation._debug?.getDragState) {
            /* no-op fallback */
          }
          if (atomPos) {
            // Derive camera world position: ArcRotateCamera may expose position directly
            const camPos = camera.position ? camera.position : (camera.getFrontPosition ? camera.getFrontPosition(0) : { x:0, y:0, z: -10 });
            const nx = atomPos.x - camPos.x;
            const ny = atomPos.y - camPos.y;
            const nz = atomPos.z - camPos.z;
            const len = Math.hypot(nx, ny, nz);
            if (len > 1e-6) {
              planeNormal = { x: nx/len, y: ny/len, z: nz/len };
              planePoint = { x: atomPos.x, y: atomPos.y, z: atomPos.z };
            }
          }
        }
      } catch { /* ignore, fallback handled below */ }
      if (!planeNormal) { planeNormal = { x:0, y:1, z:0 }; }
      if (!planePoint) { planePoint = { x:0, y:0, z:0 }; }
      const intersector = (planePoint, planeNormal) => {
        try {
          const ray = scene.createPickingRay(scene.pointerX, scene.pointerY, BABYLON.Matrix.Identity(), camera);
          const pp = new BABYLON.Vector3(planePoint.x, planePoint.y, planePoint.z);
          const pn = new BABYLON.Vector3(planeNormal.x, planeNormal.y, planeNormal.z);
          const denom = BABYLON.Vector3.Dot(pn, ray.direction);
            if (Math.abs(denom) < 1e-6) return null;
          const t = BABYLON.Vector3.Dot(pn, pp.subtract(ray.origin)) / denom;
          if (!isFinite(t) || t < 0) return null;
          const hit = ray.origin.add(ray.direction.scale(t));
          return { x: hit.x, y: hit.y, z: hit.z };
        } catch (e) {
          return null;
        }
      };
      started = beginAtomDrag(intersector, { planePoint, planeNormal });
      if (started && camera && camera.detachControl && !cameraDetachedForDrag) {
        const canvas = scene.getEngine && scene.getEngine().getRenderingCanvas ? scene.getEngine().getRenderingCanvas() : undefined;
        try { camera.detachControl(canvas); cameraDetachedForDrag = true; if (DBG) console.log('[pickingService] detach camera for drag'); } catch (e) { /* ignore */ }
        cameraLock = { alpha: camera.alpha, beta: camera.beta, radius: camera.radius, target: camera.target && { x: camera.target.x, y: camera.target.y, z: camera.target.z } };
        dragActive = true;
      }
    }
  }
  function handlePointerMove() {
    __count('pickingService#handlePointerMove');
    if (!manipulation || !cameraDetachedForDrag || !dragActive) return;
    const intersector = (planePoint, planeNormal) => {
      try {
        const ray = scene.createPickingRay(scene.pointerX, scene.pointerY, BABYLON.Matrix.Identity(), camera);
        const pp = new BABYLON.Vector3(planePoint.x, planePoint.y, planePoint.z);
        const pn = new BABYLON.Vector3(planeNormal.x, planeNormal.y, planeNormal.z);
        const denom = BABYLON.Vector3.Dot(pn, ray.direction);
        if (Math.abs(denom) < 1e-6) return null;
        const t = BABYLON.Vector3.Dot(pn, pp.subtract(ray.origin)) / denom;
        if (!isFinite(t) || t < 0) return null;
        const hit = ray.origin.add(ray.direction.scale(t));
        return { x: hit.x, y: hit.y, z: hit.z };
      } catch { return null; }
    };
  const changed = updateDragFromPointer(intersector);
    if (changed && typeof energyHook === 'function') {
      try { energyHook({ kind:'dragMove' }); } catch {}
    }
  }
  function handlePointerUp() {
    __count('pickingService#handlePointerUp');
    if (!manipulation || !cameraDetachedForDrag) return;
    endDrag();
    dragActive = false;
    if (typeof energyHook === 'function') {
      try { energyHook({ kind:'dragEnd' }); } catch {}
    }
    if (camera && camera.attachControl) {
      const canvas = scene.getEngine && scene.getEngine().getRenderingCanvas ? scene.getEngine().getRenderingCanvas() : undefined;
      try { camera.attachControl(canvas, true); if (DBG) console.log('[pickingService] reattach camera after drag'); } catch (e) { /* ignore */ }
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
    canvas && canvas.addEventListener('pointerdown', handlePointerDown);
    canvas && canvas.addEventListener('pointermove', handlePointerMove);
    canvas && canvas.addEventListener('pointerup', handlePointerUp);
    canvas && canvas.addEventListener('pointerleave', handlePointerUp);
    // Touch mapping: only attach if custom touchControls not installed
    const skipTouch = (typeof window !== 'undefined') && !!window.__MLIPVIEW_TOUCH_INSTALLED;
    if (!skipTouch) {
      const updateFromTouch = (e) => {
        try {
          const t = e.changedTouches && e.changedTouches[0];
          if (!t) return;
          const rect = typeof canvas.getBoundingClientRect === 'function' ? canvas.getBoundingClientRect() : { left:0, top:0 };
          scene.pointerX = Math.round((t.clientX ?? 0) - rect.left);
          scene.pointerY = Math.round((t.clientY ?? 0) - rect.top);
        } catch {}
      };
      const td = (e)=>{ updateFromTouch(e); try{ e.preventDefault(); e.stopPropagation(); }catch{} handlePointerDown(); };
      const tm = (e)=>{ updateFromTouch(e); try{ e.preventDefault(); e.stopPropagation(); }catch{} handlePointerMove(); };
      const tu = (e)=>{ updateFromTouch(e); try{ e.preventDefault(); e.stopPropagation(); }catch{} handlePointerUp(); };
      canvas && canvas.addEventListener('touchstart', td, { passive: false });
      canvas && canvas.addEventListener('touchmove', tm, { passive: false });
      canvas && canvas.addEventListener('touchend', tu, { passive: false });
      canvas && canvas.addEventListener('touchcancel', tu, { passive: false });
    } else if (DBG) {
      console.log('[pickingService] skipping touch handlers: custom touchControls installed');
    }
  }
  // Per-frame observer to enforce freeze while dragging
  scene.onBeforeRenderObservable && scene.onBeforeRenderObservable.add(freezeCameraFrame);
  // --- Unified interaction helpers ---
  function beginAtomDrag(intersector, opts){
    __count('pickingService#beginAtomDrag');
    if (!manipulation) return false;
    return manipulation.beginDrag(intersector, opts);
  }
  function updateDragFromPointer(intersector){
    __count('pickingService#updateDragFromPointer');
    if (!manipulation) return false;
    return manipulation.updateDrag(intersector);
  }
  function endDrag(){
    __count('pickingService#endDragUnified');
    if (!manipulation) return;
    manipulation.endDrag();
  }
  function rotateSelectedBond(deltaAngle){
    __count('pickingService#rotateSelectedBond');
    if (!manipulation) return false;
    return manipulation.rotateBond(deltaAngle);
  }
  // Expose consolidated API (supports existing pickAtPointer usage + new helpers)
  return { ...api, _debug: { get dragActive(){ return dragActive; } } };
}
