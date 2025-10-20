/** @jest-environment jsdom */

import installTouchControls from '../public/ui/touchControls.js';

function mkCanvas() {
  const c = document.createElement('canvas');
  Object.defineProperty(c, 'clientWidth', { value: 300, configurable: true });
  Object.defineProperty(c, 'clientHeight', { value: 300, configurable: true });
  c.getBoundingClientRect = () => ({ left: 0, top: 0, width: 300, height: 300 });
  return c;
}

describe('desktop rotation inputs remain intact', () => {
  test('installTouchControls does not remove camera pointer inputs on desktop', () => {
    const canvas = mkCanvas();
    const scene = { pointerX: 0, pointerY: 0 };
    const camera = {
      alpha: 1,
      beta: 1,
      radius: 10,
      attachControl: jest.fn(),
      detachControl: jest.fn(),
      inputs: {
        removeByType: jest.fn(),
        remove: jest.fn(),
        attached: { pointers: {} },
      },
    };

    // Desktop-like environment: no touch events fired. Installing should not
    // strip Babylon pointer inputs; it should just add touch listeners.
    const res = installTouchControls({ canvas, scene, camera });
    expect(res && typeof res.dispose).toBe('function');

    expect(camera.inputs.removeByType).not.toHaveBeenCalled();
    expect(camera.inputs.remove).not.toHaveBeenCalled();
  });
});
