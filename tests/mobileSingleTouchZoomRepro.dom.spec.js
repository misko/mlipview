/** @jest-environment jsdom */

import { installTouchControls } from '../public/ui/touchControls.js';

describe('mobile: single-finger should not zoom (bug repro)', () => {
  function mkCanvas() {
    const c = document.createElement('canvas');
    Object.defineProperty(c, 'clientWidth', { value: 300, configurable: true });
    Object.defineProperty(c, 'clientHeight', { value: 300, configurable: true });
    c.getBoundingClientRect = () => ({ left: 0, top: 0, width: 300, height: 300 });
    return c;
  }
  function mkScene(canvas) {
    return {
      pointerX: 0, pointerY: 0,
      onPointerObservable: { add: () => {} },
      getEngine(){ return { getRenderingCanvas(){ return canvas; } }; },
      onBeforeRenderObservable: { add: () => {} }
    };
  }

  test('without protective detach, single touch causes radius to change (repro), with fix it stays constant', () => {
    const canvas = mkCanvas();
    document.body.appendChild(canvas);

    // Camera stub that zooms on pointermove if attached (simulates Babylon input still active)
    let attached = false;
    const onMove = () => { if (attached) camera.radius += 1; };
    const camera = {
      alpha: 1, beta: 1, radius: 10,
      inputs: { attached: { pointers: {} }, removeByType: jest.fn(), remove: jest.fn() },
      attachControl: jest.fn(() => { canvas.addEventListener('pointermove', onMove); attached = true; }),
      detachControl: jest.fn(() => { canvas.removeEventListener('pointermove', onMove); attached = false; })
    };
    // Attach once like the app does
    camera.attachControl(canvas, true);

    const scene = mkScene(canvas);

    // Install touch controls - prior to fix this won't stop the camera from zooming on single finger
    installTouchControls({ canvas, scene, camera, picking: { _debug: { get dragActive(){ return false; } } } });

    const start = new Event('touchstart', { bubbles: true, cancelable: true });
    start.changedTouches = [{ identifier: 0, clientX: 100, clientY: 100 }];
    canvas.dispatchEvent(start);

    const move = new Event('touchmove', { bubbles: true, cancelable: true });
    move.changedTouches = [{ identifier: 0, clientX: 140, clientY: 120 }];
    canvas.dispatchEvent(move);

    // If the underlying camera pointers input is still active, our synthetic pointermove will cause zoom (radius++).
    // The fix ensures we temporarily detach camera during touch, so radius must remain the same.
    expect(camera.radius).toBe(10);

    const end = new Event('touchend', { bubbles: true, cancelable: true });
    end.changedTouches = [{ identifier: 0, clientX: 140, clientY: 120 }];
    canvas.dispatchEvent(end);
  });
});
