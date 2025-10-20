/** @jest-environment jsdom */

import { installTouchControls } from '../public/ui/touchControls.js';

describe('mobile: camera detach/reattach during touch gestures', () => {
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
      getEngine() {
        return {
          getRenderingCanvas() {
            return canvas;
          },
        };
      },
    };
  }

  function touchEvent(type, points) {
    const ev = new Event(type, { bubbles: true, cancelable: true });
    ev.changedTouches = points.map((p, i) => ({ identifier: i, clientX: p.x, clientY: p.y }));
    return ev;
  }

  test('single-finger drag detaches camera, rotates alpha, preserves radius, then reattaches', () => {
    const canvas = mkCanvas();
    document.body.appendChild(canvas);

    // Simulate an attached camera that would normally zoom on pointermove if left attached
    let isAttached = false;
    const onPointerMove = () => {
      if (isAttached) camera.radius += 1000;
    }; // extreme zoom to catch any leakage

    const camera = {
      alpha: 1,
      beta: 1,
      radius: 10,
      inertialRadiusOffset: 0,
      inputs: { attached: { pointers: {} }, removeByType: jest.fn(), remove: jest.fn() },
      attachControl: jest.fn(() => {
        isAttached = true;
        canvas.addEventListener('pointermove', onPointerMove);
      }),
      detachControl: jest.fn(() => {
        isAttached = false;
        canvas.removeEventListener('pointermove', onPointerMove);
      }),
    };

    // App would usually attach once at init
    camera.attachControl(canvas, true);

    const scene = mkScene(canvas);
    const picking = {
      _debug: {
        get dragActive() {
          return false;
        },
      },
    };

    installTouchControls({ canvas, scene, camera, picking });

    // Begin single-finger gesture
    canvas.dispatchEvent(touchEvent('touchstart', [{ x: 100, y: 100 }]));

    // Should have detached camera controls to avoid default zoom behavior
    expect(camera.detachControl).toHaveBeenCalled();

    const a0 = camera.alpha;
    const r0 = camera.radius;

    // Move finger - should rotate (alpha != a0) and not zoom (radius === r0)
    canvas.dispatchEvent(touchEvent('touchmove', [{ x: 140, y: 120 }]));

    expect(camera.alpha).not.toBe(a0);
    expect(camera.radius).toBe(r0);

    // End gesture
    canvas.dispatchEvent(touchEvent('touchend', [{ x: 140, y: 120 }]));

    // After touch ends, controls should be reattached (no drag active)
    expect(camera.attachControl).toHaveBeenCalled();
  });
});
