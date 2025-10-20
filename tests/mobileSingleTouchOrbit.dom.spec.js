/** @jest-environment jsdom */

import { installTouchControls } from '../public/ui/touchControls.js';

function mkCanvas() {
  const c = document.createElement('canvas');
  Object.defineProperty(c, 'clientWidth', { value: 300, configurable: true });
  Object.defineProperty(c, 'clientHeight', { value: 300, configurable: true });
  c.getBoundingClientRect = () => ({ left: 0, top: 0, width: 300, height: 300 });
  return c;
}

function mkScene(canvas) {
  return {
    pointerX: 0,
    pointerY: 0,
    onPointerObservable: { add: () => {} },
    getEngine() {
      return {
        getRenderingCanvas() {
          return canvas;
        },
      };
    },
    onBeforeRenderObservable: { add: () => {} },
  };
}

describe('mobile: single-finger orbit rotation', () => {
  test('single touch drag rotates camera alpha and preserves radius', () => {
    const canvas = mkCanvas();
    document.body.appendChild(canvas);

    const camera = { alpha: 1, beta: 1, radius: 10 };
    const scene = mkScene(canvas);

    // Picking not needed for empty-space orbit; ensure dragActive is false if referenced
    const picking = {
      _debug: {
        get dragActive() {
          return false;
        },
      },
    };

    installTouchControls({ canvas, scene, camera, picking });

    const start = new Event('touchstart', { bubbles: true, cancelable: true });
    start.changedTouches = [{ identifier: 0, clientX: 100, clientY: 100 }];
    canvas.dispatchEvent(start);

    const move = new Event('touchmove', { bubbles: true, cancelable: true });
    move.changedTouches = [{ identifier: 0, clientX: 140, clientY: 120 }];
    canvas.dispatchEvent(move);

    // Alpha should have changed due to dx > 0; radius should remain the same for single-finger rotate
    expect(camera.alpha).not.toBe(1);
    expect(camera.radius).toBe(10);

    const end = new Event('touchend', { bubbles: true, cancelable: true });
    end.changedTouches = [{ identifier: 0, clientX: 140, clientY: 120 }];
    canvas.dispatchEvent(end);
  });
});
