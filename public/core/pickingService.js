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
export function createPickingService(
  scene,
  view,
  selectionService,
  { manipulation, camera, energyHook, bondScrollStep } = {}
) {
  __count('pickingService#createPickingService');
  let cameraDetachedForDrag = false;
  let cameraLock = null; // { alpha,beta,radius,target }
  let dragActive = false;
  const DBG = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_TOUCH;
  const PDBG = typeof window !== 'undefined' && !!window.__MLIPVIEW_DEBUG_PICK;
  let interactionsEnabled = true;
  // De-dupe guards
  // 1) Exact event de-dup across capture/bubble/observable paths
  const __handledPointerDownEvents = typeof WeakSet !== 'undefined' ? new WeakSet() : null;
  // 2) Spatial-temporal guard as a safety net
  // Last processed pointerdown for spatial-temporal de-duplication; null until first real event
  let __lastDown = null;
  // Public interaction surface (populated below) returned at end
  const api = {
    pickAtPointer,
    beginAtomDrag,
    updateDragFromPointer,
    endDrag,
    rotateSelectedBond,
    freezeCameraFrame,
    selection: selectionService.get,
    setEnabled(on = true) {
      const next = !!on;
      if (!next && dragActive) {
        const prev = interactionsEnabled;
        interactionsEnabled = true;
        try { handlePointerUp(); } catch { }
        interactionsEnabled = prev;
      }
      interactionsEnabled = next;
      return interactionsEnabled;
    },
    isEnabled() {
      return interactionsEnabled;
    },
  };
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
    if (!interactionsEnabled) return null;
    if (PDBG) {
      try {
        console.log('[pick] at', {
          x: scene.pointerX,
          y: scene.pointerY,
          hasMulti: typeof scene.multiPick === 'function',
        });
      } catch {}
    }
    const pick = scene.pick(scene.pointerX, scene.pointerY);
    if (!pick || !pick.hit) return null;
    // Resolve atom first, then bond
    const atom = view.resolveAtomPick(pick);
    if (atom) {
      if (PDBG)
        try {
          console.log('[pick] singlePick atom', atom);
        } catch {}
      return atom;
    }
    const bond = view.resolveBondPick(pick);
    if (bond) {
      if (PDBG)
        try {
          console.log('[pick] singlePick bond', bond);
        } catch {}
      return bond;
    }
    return null;
  }
  function handlePointerDown(e) {
    __count('pickingService#handlePointerDown');
    if (!interactionsEnabled) return;
    // Hard de-dup: if we've already processed this DOM event, skip.
    try {
      if (e && __handledPointerDownEvents) {
        if (__handledPointerDownEvents.has(e)) return;
        __handledPointerDownEvents.add(e);
      }
    } catch {}
    try {
      const now =
        typeof performance !== 'undefined' && performance.now ? performance.now() : Date.now();
      if (__lastDown) {
        const dx = Math.abs((scene.pointerX || 0) - (__lastDown.x || 0));
        const dy = Math.abs((scene.pointerY || 0) - (__lastDown.y || 0));
        if (dx <= 1 && dy <= 1 && now - (__lastDown.t || 0) < 150) {
          if (PDBG) console.log('[pick] dedupe pointerdown');
          return;
        }
      }
      __lastDown = { x: scene.pointerX || 0, y: scene.pointerY || 0, t: now };
    } catch {}
    const res = pickAtPointer();
    const isRightButton = e && typeof e.button === 'number' && e.button === 2;
    if (!res) {
      if (PDBG)
        try {
          console.log('[pick] no hit; clearing selection');
        } catch {}
      if (!isRightButton && selectionService.get().kind !== 'atom' && view && view._debugAutoSelectFirstOnEmpty) {
        selectionService.clickAtom(0);
      } else {
        selectionService.clear();
        return;
      }
    }
    if (isRightButton) {
      if (res && res.kind === 'atom') {
        try {
          if (typeof e?.preventDefault === 'function') e.preventDefault();
          if (typeof e?.stopPropagation === 'function') e.stopPropagation();
          if (typeof e?.stopImmediatePropagation === 'function') e.stopImmediatePropagation();
        } catch {}
        try {
          const api = typeof window !== 'undefined' ? window.viewerApi || null : null;
          if (api && typeof api.removeAtomByIndex === 'function') {
            Promise.resolve(api.removeAtomByIndex(res.index)).catch(() => {});
          }
        } catch {}
      }
      return;
    }
    if (res.kind === 'atom') {
      if (PDBG)
        try {
          console.log('[select] atom', res.index);
        } catch {}
      selectionService.clickAtom(res.index);
    } else if (res.kind === 'bond') {
      if (PDBG)
        try {
          console.log('[select] bond', res);
        } catch {}
      selectionService.clickBond(res);
    }
    // After selection, if atom selected and manipulation present, initiate drag baseline
    if (manipulation && selectionService.get().kind === 'atom') {
      const atomSel = selectionService.get();
      let started = false;
      // Compute camera-aligned drag plane: normal = (atomPos - cameraPos)
      let planePoint = null;
      let planeNormal = null;
      try {
        if (camera && atomSel && atomSel.kind === 'atom') {
          const atomPos = manipulation.molState
            ? manipulation.molState.positions[atomSel.data.index]
            : view && view._internals && view._internals.highlight && view._internals.highlight.atom
              ? view._internals.highlight.atom.position
              : null;
          // Fallback: take position from molState via selection if available
          if (!atomPos && manipulation && manipulation._debug?.getDragState) {
            /* no-op fallback */
          }
          if (atomPos) {
            // Derive camera world position: ArcRotateCamera may expose position directly
            const camPos = camera.position
              ? camera.position
              : camera.getFrontPosition
                ? camera.getFrontPosition(0)
                : { x: 0, y: 0, z: -10 };
            const nx = atomPos.x - camPos.x;
            const ny = atomPos.y - camPos.y;
            const nz = atomPos.z - camPos.z;
            const len = Math.hypot(nx, ny, nz);
            if (len > 1e-6) {
              planeNormal = { x: nx / len, y: ny / len, z: nz / len };
              planePoint = { x: atomPos.x, y: atomPos.y, z: atomPos.z };
            }
          }
        }
      } catch {
        /* ignore, fallback handled below */
      }
      if (!planeNormal) {
        planeNormal = { x: 0, y: 1, z: 0 };
      }
      if (!planePoint) {
        planePoint = { x: 0, y: 0, z: 0 };
      }
      const intersector = (planePoint, planeNormal) => {
        try {
          const ray = scene.createPickingRay(
            scene.pointerX,
            scene.pointerY,
            BABYLON.Matrix.Identity(),
            camera
          );
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
        const canvas =
          scene.getEngine && scene.getEngine().getRenderingCanvas
            ? scene.getEngine().getRenderingCanvas()
            : undefined;
        try {
          camera.detachControl(canvas);
          cameraDetachedForDrag = true;
          if (DBG) console.log('[pickingService] detach camera for drag');
        } catch (e) {
          /* ignore */
        }
        cameraLock = {
          alpha: camera.alpha,
          beta: camera.beta,
          radius: camera.radius,
          target: camera.target && { x: camera.target.x, y: camera.target.y, z: camera.target.z },
        };
        dragActive = true;
      }
    }
  }
  function handlePointerMove() {
    __count('pickingService#handlePointerMove');
    if (!interactionsEnabled) return;
    if (!manipulation || !cameraDetachedForDrag || !dragActive) return;
    const intersector = (planePoint, planeNormal) => {
      try {
        const ray = scene.createPickingRay(
          scene.pointerX,
          scene.pointerY,
          BABYLON.Matrix.Identity(),
          camera
        );
        const pp = new BABYLON.Vector3(planePoint.x, planePoint.y, planePoint.z);
        const pn = new BABYLON.Vector3(planeNormal.x, planeNormal.y, planeNormal.z);
        const denom = BABYLON.Vector3.Dot(pn, ray.direction);
        if (Math.abs(denom) < 1e-6) return null;
        const t = BABYLON.Vector3.Dot(pn, pp.subtract(ray.origin)) / denom;
        if (!isFinite(t) || t < 0) return null;
        const hit = ray.origin.add(ray.direction.scale(t));
        return { x: hit.x, y: hit.y, z: hit.z };
      } catch {
        return null;
      }
    };
    const changed = updateDragFromPointer(intersector);
    if (changed && typeof energyHook === 'function') {
      try {
        energyHook({ kind: 'dragMove' });
      } catch {}
    }
  }
  function handlePointerUp() {
    __count('pickingService#handlePointerUp');
    if (!interactionsEnabled) return;
    if (!manipulation || !cameraDetachedForDrag) return;
    endDrag();
    dragActive = false;
    if (typeof energyHook === 'function') {
      try {
        energyHook({ kind: 'dragEnd' });
      } catch {}
    }
    if (camera && camera.attachControl) {
      const canvas =
        scene.getEngine && scene.getEngine().getRenderingCanvas
          ? scene.getEngine().getRenderingCanvas()
          : undefined;
      try {
        camera.attachControl(canvas, true);
        if (DBG) console.log('[pickingService] reattach camera after drag');
      } catch (e) {
        /* ignore */
      }
    }
    cameraDetachedForDrag = false;
    cameraLock = null;
  }
  // Prefer DOM listeners when a rendering canvas exists; fall back to Babylon pointer observable otherwise.
  let usedDomListeners = false;
  if (scene.getEngine && scene.getEngine().getRenderingCanvas) {
    const canvas = scene.getEngine().getRenderingCanvas();
    const effectiveBondScrollStep =
      typeof bondScrollStep === 'number' && Number.isFinite(bondScrollStep) && bondScrollStep > 0
        ? bondScrollStep
        : Math.PI / 36;
    // Avoid duplicate pointerdown handling: scene.onPointerObservable already listens for POINTERDOWN above.
    // We still attach DOM listeners to support environments where Babylon doesn't pipe DOM -> onPointerObservable (tests/jsdom),
    // and to track pointer coordinates for drag.
    const updateFromPointer = (e) => {
      try {
        // Skip synthetic pointer events emitted by touchControls; that layer already
        // set scene.pointerX/Y from actual touch coordinates just before dispatch.
        if (e && e.__mlip_synthetic_from_touch) return;
        const rect =
          typeof canvas.getBoundingClientRect === 'function'
            ? canvas.getBoundingClientRect()
            : { left: 0, top: 0 };
        const cx = e && typeof e.clientX === 'number' ? e.clientX : 0;
        const cy = e && typeof e.clientY === 'number' ? e.clientY : 0;
        scene.pointerX = Math.round(cx - rect.left);
        scene.pointerY = Math.round(cy - rect.top);
      } catch {}
    };
    const pd = (e) => {
      try {
        e.__mlip_handled_byPickingService = true;
      } catch {}
      updateFromPointer(e);
      handlePointerDown(e);
    };
    const pm = (e) => {
      updateFromPointer(e);
      handlePointerMove();
    };
    const pu = (e) => {
      updateFromPointer(e);
      handlePointerUp();
    };
    canvas && canvas.addEventListener('pointerdown', pd);
    canvas && canvas.addEventListener('pointermove', pm);
    canvas && canvas.addEventListener('pointerup', pu);
    canvas && canvas.addEventListener('pointerleave', pu);
    canvas && canvas.addEventListener('contextmenu', (e) => {
      if (typeof e?.preventDefault === 'function') e.preventDefault();
      if (typeof e?.stopPropagation === 'function') e.stopPropagation();
      if (typeof e?.stopImmediatePropagation === 'function') e.stopImmediatePropagation();
      return false;
    });
    const onWheel = (e) => {
      try {
        if (!manipulation || typeof manipulation.rotateBond !== 'function') return;
        const sel = selectionService?.get?.();
        if (!sel || sel.kind !== 'bond') return;
        if (typeof e?.preventDefault === 'function') e.preventDefault();
        if (typeof e?.stopPropagation === 'function') e.stopPropagation();
        if (typeof e?.stopImmediatePropagation === 'function') e.stopImmediatePropagation();
        const deltaY = typeof e?.deltaY === 'number' ? e.deltaY : 0;
        if (!Number.isFinite(deltaY) || deltaY === 0) return;
        const scale = -deltaY / 120;
        const angle = scale * effectiveBondScrollStep;
        if (!Number.isFinite(angle) || Math.abs(angle) < 1e-5) return;
        const rotated = manipulation.rotateBond(angle);
        if (rotated && typeof energyHook === 'function') {
          try {
            energyHook({ kind: 'bondRotate' });
          } catch {}
        }
      } catch {}
    };
    canvas && canvas.addEventListener('wheel', onWheel, { passive: false, capture: true });
    usedDomListeners = !!canvas;
    // Removed document-level capture fallback to reduce complexity; canvas listeners suffice in browsers/tests.
    // Touch mapping: only attach if custom touchControls not installed
    const skipTouch =
      typeof window !== 'undefined' &&
      (!!window.__MLIPVIEW_TOUCH_INSTALLED || !!window.__MLIPVIEW_NO_TOUCH);
    if (!skipTouch) {
      if (PDBG)
        try {
          console.log('[pickingService] attaching built-in touch listeners');
        } catch {}
      const updateFromTouch = (e) => {
        try {
          const t = e.changedTouches && e.changedTouches[0];
          if (!t) return;
          const rect =
            typeof canvas.getBoundingClientRect === 'function'
              ? canvas.getBoundingClientRect()
              : { left: 0, top: 0 };
          scene.pointerX = Math.round((t.clientX ?? 0) - rect.left);
          scene.pointerY = Math.round((t.clientY ?? 0) - rect.top);
        } catch {}
      };
      const td = (e) => {
        updateFromTouch(e);
        try {
          e.preventDefault();
          e.stopPropagation();
        } catch {}
        handlePointerDown();
      };
      const tm = (e) => {
        updateFromTouch(e);
        try {
          e.preventDefault();
          e.stopPropagation();
        } catch {}
        handlePointerMove();
      };
      const tu = (e) => {
        updateFromTouch(e);
        try {
          e.preventDefault();
          e.stopPropagation();
        } catch {}
        handlePointerUp();
      };
      canvas && canvas.addEventListener('touchstart', td, { passive: false });
      canvas && canvas.addEventListener('touchmove', tm, { passive: false });
      canvas && canvas.addEventListener('touchend', tu, { passive: false });
      canvas && canvas.addEventListener('touchcancel', tu, { passive: false });
    } else if (PDBG || DBG) {
      console.log(
        '[pickingService] skipping touch handlers:',
        window.__MLIPVIEW_NO_TOUCH ? 'noTouch mode' : 'custom touchControls installed'
      );
    }
  }
  // Always register Babylon pointer observable handler as a secondary path so unit tests
  // can trigger POINTERDOWN via scene.onPointerObservable even if DOM listeners are active.
  scene.onPointerObservable.add((pi) => {
    if (pi.type === BABYLON.PointerEventTypes.POINTERDOWN) {
      // If this POINTERDOWN originated from the same DOM event our own listener already handled,
      // skip to avoid double-processing (which can cycle/clear bond selection).
      try {
        if (pi.event && pi.event.__mlip_handled_byPickingService) return;
      } catch {}
      handlePointerDown(pi.event);
    }
  });
  // Per-frame observer to enforce freeze while dragging
  scene.onBeforeRenderObservable && scene.onBeforeRenderObservable.add(freezeCameraFrame);
  // --- Unified interaction helpers ---
  function beginAtomDrag(intersector, opts) {
    __count('pickingService#beginAtomDrag');
    if (!manipulation) return false;
    return manipulation.beginDrag(intersector, opts);
  }
  function updateDragFromPointer(intersector) {
    __count('pickingService#updateDragFromPointer');
    if (!manipulation) return false;
    return manipulation.updateDrag(intersector);
  }
  function endDrag() {
    __count('pickingService#endDragUnified');
    if (!manipulation) return;
    manipulation.endDrag();
  }
  function rotateSelectedBond(deltaAngle) {
    __count('pickingService#rotateSelectedBond');
    if (!manipulation) return false;
    return manipulation.rotateBond(deltaAngle);
  }
  // Expose consolidated API (supports existing pickAtPointer usage + new helpers)
  return {
    ...api,
    _debug: {
      get dragActive() {
        return dragActive;
      },
    },
  };
}
