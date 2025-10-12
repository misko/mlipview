/** @jest-environment jsdom */

import { installTouchControls } from '../public/ui/touchControls.js';

describe('mobile gestures', () => {
	function mkCanvas() {
		const c = document.createElement('canvas');
		Object.defineProperty(c, 'clientWidth', { value: 200, configurable: true });
		Object.defineProperty(c, 'clientHeight', { value: 200, configurable: true });
		c.getBoundingClientRect = () => ({ left: 0, top: 0, width: 200, height: 200 });
		return c;
	}
	function mkScene(canvas) {
		return {
			pointerX: 0,
			pointerY: 0,
			onPointerObservable: { add: jest.fn() },
			onBeforeRenderObservable: { add: jest.fn() },
			getEngine() { return { getRenderingCanvas() { return canvas; } }; },
			pick: jest.fn(() => ({ hit: false })),
			createPickingRay: jest.fn(() => ({ origin: { x: 0, y: 0, z: -5 }, direction: { x: 0, y: 0, z: 1 } }))
		};
	}
	function touchEvent(type, points) {
		const ev = new Event(type, { bubbles: true, cancelable: true });
		ev.changedTouches = points.map((p, i) => ({ identifier: i, clientX: p.x, clientY: p.y }));
		return ev;
	}

	test('pinch zoom updates camera radius', () => {
		const canvas = mkCanvas(); document.body.appendChild(canvas);
		const camera = { alpha: 1, beta: 1, radius: 10, inertialRadiusOffset: 0, position: { x: 0, y: 0, z: -10 } };
		const scene = mkScene(canvas);
		const picking = { _debug: { dragActive: false } };
		installTouchControls({ canvas, scene, camera, picking });

		const r0 = camera.radius;
		// start with two fingers
		canvas.dispatchEvent(touchEvent('touchstart', [{ x: 80, y: 100 }, { x: 120, y: 100 }]));
		// move fingers apart
		canvas.dispatchEvent(touchEvent('touchmove', [{ x: 40, y: 100 }, { x: 160, y: 100 }]));
		expect(camera.radius).toBeLessThan(r0);
	});

	test('single-finger drag orbits camera (alpha changes)', () => {
		const canvas = mkCanvas(); document.body.appendChild(canvas);
		const camera = { alpha: 1, beta: 1, radius: 10, inertialRadiusOffset: 0, position: { x: 0, y: 0, z: -10 } };
		const scene = mkScene(canvas);
		const picking = { _debug: { dragActive: false } };
		installTouchControls({ canvas, scene, camera, picking });

		canvas.dispatchEvent(touchEvent('touchstart', [{ x: 100, y: 100 }]));
		const a0 = camera.alpha;
		canvas.dispatchEvent(touchEvent('touchmove', [{ x: 150, y: 100 }]));
		expect(camera.alpha).not.toBe(a0);
	});
});
