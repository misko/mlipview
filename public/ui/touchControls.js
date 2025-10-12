// Touch controls layer
// - Maps touch gestures to camera orbit (single-finger) and zoom (pinch)
// - Updates scene.pointerX/Y and synthesizes pointer events so picking/drag logic stays intact

export function installTouchControls({ canvas, scene, camera, picking } = {}) {
	if (!canvas || !scene) return;
	try { if (typeof window !== 'undefined' && window.__MLIPVIEW_NO_TOUCH) { console.warn('[touchControls] install aborted by __MLIPVIEW_NO_TOUCH'); return; } } catch {}
	// Mark globally so other systems (like pickingService) can avoid double-attaching touch listeners
	try { if (typeof window !== 'undefined') window.__MLIPVIEW_TOUCH_INSTALLED = true; } catch {}
	const DBG = (typeof window !== 'undefined') && !!window.__MLIPVIEW_DEBUG_TOUCH;
	const rectOf = () => (typeof canvas.getBoundingClientRect === 'function' ? canvas.getBoundingClientRect() : { left: 0, top: 0 });

  // Prevent native scrolling/zoom from interfering with gestures
	try { if (canvas && canvas.style) { canvas.style.touchAction = 'none'; canvas.style.webkitUserSelect = 'none'; } } catch {}

  // Gesture tracking state
	const touches = new Map(); // id -> {x,y}
	let mode = 'none'; // 'none' | 'rotate' | 'pinch'
	let lastSingle = null; // {x,y}
	let pinchRef = null; // {d0, r0}

	let cameraDetached = false;
	const getDragging = () => !!(picking && picking._debug && picking._debug.dragActive);
  const detachCamera = () => {
		try {
			if (camera && camera.detachControl && !cameraDetached) {
				camera.detachControl(canvas);
				cameraDetached = true;
				if (DBG) console.log('[touchControls] detachCamera on touchstart');
			}
		} catch {}
	};
  const attachCameraIfIdle = () => {
		try {
			if (camera && camera.attachControl && cameraDetached && !getDragging()) {
				camera.attachControl(canvas, true);
				cameraDetached = false;
				if (DBG) console.log('[touchControls] attachCamera after touchend');
			}
		} catch {}
	};

  const updatePointer = (t) => {
		try {
			const r = rectOf();
			scene.pointerX = Math.round((t.clientX ?? 0) - r.left);
			scene.pointerY = Math.round((t.clientY ?? 0) - r.top);
		} catch {}
	};
  const synth = (type) => {
		try {
			const ev = new Event(type, { bubbles: true, cancelable: true });
			try { ev.__mlip_synthetic_from_touch = true; } catch {}
			canvas.dispatchEvent(ev);
		} catch {
			// Fallback without options for older environments
			try { const ev = new Event(type); try { ev.__mlip_synthetic_from_touch = true; } catch {}; canvas.dispatchEvent(ev); } catch {}
		}
	};

  const onStart = (e) => {
		if (!e || !e.changedTouches || e.changedTouches.length === 0) return;
		// Track touches
		for (const t of e.changedTouches) {
			touches.set(t.identifier, { x: t.clientX, y: t.clientY });
		}
		const t = e.changedTouches[0];
		updatePointer(t);
		// Determine mode
		if (touches.size >= 2) {
			// Pinch start
			const arr = Array.from(touches.values());
			const dx = arr[1].x - arr[0].x; const dy = arr[1].y - arr[0].y;
			const d0 = Math.hypot(dx, dy) || 1;
			pinchRef = { d0, r0: camera && typeof camera.radius === 'number' ? camera.radius : 10 };
			mode = 'pinch';
			if (DBG) console.log('[touchControls] mode=pinch start', { d0, r0: pinchRef.r0 });
		} else if (touches.size === 1) {
			// Rotate start
			lastSingle = { x: t.clientX, y: t.clientY };
			mode = 'rotate';
			if (DBG) console.log('[touchControls] mode=rotate start at', lastSingle);
		}
		// While touch is active, detach camera controls to avoid native zoom/pan from pointers input
		detachCamera();
		// Prevent default scrolling/zoom
		try { e.preventDefault(); e.stopPropagation(); } catch {}
		// Kick picking via pointerdown after pointer coords are updated so picking reads correct scene.pointerX/Y
		synth('pointerdown');
	};
  const onMove = (e) => {
		if (!e || !e.changedTouches || e.changedTouches.length === 0) return;
		// Update tracked touches
		for (const t of e.changedTouches) {
			touches.set(t.identifier, { x: t.clientX, y: t.clientY });
		}
		const t = e.changedTouches[0];
		updatePointer(t);
		// Apply gesture effects
		if (mode === 'pinch' && touches.size >= 2) {
			const arr = Array.from(touches.values());
			const dx = arr[1].x - arr[0].x; const dy = arr[1].y - arr[0].y;
			const d = Math.hypot(dx, dy) || 1;
			if (pinchRef && camera) {
				// Increase distance -> zoom in (smaller radius)
				const ratio = Math.max(0.05, Math.min(20, pinchRef.d0 / d));
				const next = pinchRef.r0 * ratio;
				if (typeof camera.inertialRadiusOffset === 'number') camera.inertialRadiusOffset = 0;
				camera.radius = next;
				if (DBG) console.log('[touchControls] pinch move', { ratio, radius: camera.radius });
			}
		} else if (mode === 'rotate' && touches.size === 1) {
			// Only rotate if not actively dragging an atom
			const isDragging = !!(picking && picking._debug && picking._debug.dragActive);
			if (!isDragging && camera && lastSingle) {
				const dx = t.clientX - lastSingle.x;
				const dy = t.clientY - lastSingle.y;
				// Tune sensitivity lightly
				const sensitivity = 0.01;
				if (typeof camera.inertialRadiusOffset === 'number') camera.inertialRadiusOffset = 0; // suppress zoom
				if (typeof camera.alpha === 'number') camera.alpha += dx * sensitivity;
				if (typeof camera.beta === 'number') camera.beta += dy * sensitivity * 0.5;
				if (DBG) console.log('[touchControls] rotate move', { dx, dy, alpha: camera.alpha, beta: camera.beta });
				lastSingle = { x: t.clientX, y: t.clientY };
			}
		}
		try { e.preventDefault(); e.stopPropagation(); } catch {}
		synth('pointermove');
	};
  const onEnd = (e) => {
		if (e && e.changedTouches) {
			for (const t of e.changedTouches) touches.delete(t.identifier);
			if (e.changedTouches[0]) updatePointer(e.changedTouches[0]);
		}
		// Reset modes when touches drop below thresholds
		if (touches.size < 2) pinchRef = null;
		if (touches.size === 0) { if (DBG) console.log('[touchControls] gesture end'); mode = 'none'; lastSingle = null; }
		try { e.preventDefault(); e.stopPropagation(); } catch {}
		synth('pointerup');
		// Reattach camera controls when gesture fully ends and no drag is active
		if (touches.size === 0) attachCameraIfIdle();
	};

	canvas.addEventListener('touchstart', onStart, { passive: false });
	canvas.addEventListener('touchmove', onMove, { passive: false });
	canvas.addEventListener('touchend', onEnd, { passive: false });
	canvas.addEventListener('touchcancel', onEnd, { passive: false });

	// Do NOT strip Babylon pointer inputs on desktop; it breaks mouse rotation.
	// On touch devices, we rely on temporary camera detach/attach during gestures.
	try {
		const isBrowser = (typeof window !== 'undefined');
		const touchCapable = isBrowser && (('ontouchstart' in window) || (typeof navigator !== 'undefined' && navigator.maxTouchPoints > 0));
		if (touchCapable && camera && camera.inputs) {
			// Historically we attempted to remove pointer inputs here, but detaching the camera
			// during touch is sufficient. We keep this as a no-op to avoid side-effects.
			// If a future device misbehaves, enable targeted removals behind a debug flag.
			if (isBrowser && window.__MLIPVIEW_FORCE_REMOVE_POINTER_INPUTS === true) {
				if (typeof camera.inputs.removeByType === 'function') {
					try { camera.inputs.removeByType('ArcRotateCameraPointersInput'); } catch {}
					try { camera.inputs.removeByType('FreeCameraPointersInput'); } catch {}
					try { camera.inputs.removeByType('pointers'); } catch {}
				}
				if (typeof camera.inputs.remove === 'function') {
					try { camera.inputs.remove(camera.inputs.attached && camera.inputs.attached.pointers); } catch {}
				}
			}
		}
	} catch {}

	return {
		dispose() {
			canvas.removeEventListener('touchstart', onStart);
			canvas.removeEventListener('touchmove', onMove);
			canvas.removeEventListener('touchend', onEnd);
			canvas.removeEventListener('touchcancel', onEnd);
		}
	};
}

// Provide a default export for interop with different module systems (Jest/Babel)
export default installTouchControls;

