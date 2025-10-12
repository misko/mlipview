// Minimal touch controls layer
// Purpose: bridge touch events to pointer events so existing picking/drag logic works.
// We keep this lightweight: update scene.pointerX/Y and synthesize pointerdown/move/up
// on the canvas, which the pickingService listens to.

export function installTouchControls({ canvas, scene, camera, picking } = {}) {
	if (!canvas || !scene) return;
	const rectOf = () => (typeof canvas.getBoundingClientRect === 'function' ? canvas.getBoundingClientRect() : { left: 0, top: 0 });

	// Gesture tracking state
	const touches = new Map(); // id -> {x,y}
	let mode = 'none'; // 'none' | 'rotate' | 'pinch'
	let lastSingle = null; // {x,y}
	let pinchRef = null; // {d0, r0}

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
			canvas.dispatchEvent(ev);
		} catch {
			// Fallback without options for older environments
			try { canvas.dispatchEvent(new Event(type)); } catch {}
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
		} else if (touches.size === 1) {
			// Rotate start
			lastSingle = { x: t.clientX, y: t.clientY };
			mode = 'rotate';
		}
		// Prevent default scrolling/zoom
		try { e.preventDefault(); e.stopPropagation(); } catch {}
		// Kick picking via pointerdown
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
		if (touches.size === 0) { mode = 'none'; lastSingle = null; }
		try { e.preventDefault(); e.stopPropagation(); } catch {}
		synth('pointerup');
	};

	canvas.addEventListener('touchstart', onStart, { passive: false });
	canvas.addEventListener('touchmove', onMove, { passive: false });
	canvas.addEventListener('touchend', onEnd, { passive: false });
	canvas.addEventListener('touchcancel', onEnd, { passive: false });

	// Detach default camera pointer inputs during touch drags if available (optional safe guard)
	try {
			if (camera && camera.inputs) {
				// Prefer removeByType if available so tests can assert it
				if (typeof camera.inputs.removeByType === 'function') {
					try { camera.inputs.removeByType('ArcRotateCameraPointersInput'); } catch {}
					try { camera.inputs.removeByType('FreeCameraPointersInput'); } catch {}
					try { camera.inputs.removeByType('pointers'); } catch {}
				}
				if (typeof camera.inputs.remove === 'function') {
					// Also attempt direct removal of attached pointer input
					try { camera.inputs.remove(camera.inputs.attached && camera.inputs.attached.pointers); } catch {}
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

