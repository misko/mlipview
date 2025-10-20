/** @jest-environment jsdom */

import { installTouchControls } from '../public/ui/touchControls.js';

describe('mobile: custom touch controls suppress default pointer zoom', () => {
  function mkCanvas() {
    const c = document.createElement('canvas');
    c.getBoundingClientRect = () => ({ left: 0, top: 0 });
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
      pick: jest.fn(() => ({ hit: false })),
      createPickingRay: jest.fn(() => ({
        origin: { x: 0, y: 0, z: -5 },
        direction: { x: 0, y: 0, z: 1 },
      })),
    };
  }
  function touch(type, x, y) {
    const ev = new Event(type, { bubbles: true, cancelable: true });
    ev.changedTouches = [{ identifier: 0, clientX: x, clientY: y }];
    return ev;
  }

  test('single finger rotate does not change radius even if default inputs exist', () => {
    const canvas = mkCanvas();
    document.body.appendChild(canvas);
    const camera = {
      alpha: 1,
      beta: 1,
      radius: 10,
      inertialRadiusOffset: 0,
      attachControl: jest.fn(),
      detachControl: jest.fn(),
      inputs: { attached: { pointers: {} }, remove: jest.fn(), removeByType: jest.fn() },
    };
    const scene = mkScene(canvas);
    const picking = { _debug: { dragActive: false } };
    installTouchControls({ canvas, scene, camera, picking });

    canvas.dispatchEvent(touch('touchstart', 100, 100));
    const a0 = camera.alpha;
    const r0 = camera.radius;
    canvas.dispatchEvent(touch('touchmove', 150, 100));
    canvas.dispatchEvent(touch('touchend', 150, 100));

    expect(camera.alpha).not.toBe(a0);
    expect(camera.radius).toBe(r0);
    // New behavior: we detach Babylon camera controls during touch and re-attach after.
    expect(camera.detachControl).toHaveBeenCalled();
    expect(camera.attachControl).toHaveBeenCalled();
  });
});
