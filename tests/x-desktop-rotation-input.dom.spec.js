/** @jest-environment jsdom */

import installTouchControls from '../public/ui/touchControls.js';

function makeCanvas() {
  const canvas = document.createElement('canvas');
  Object.defineProperty(canvas, 'clientWidth', { value: 320, configurable: true });
  Object.defineProperty(canvas, 'clientHeight', { value: 240, configurable: true });
  canvas.getBoundingClientRect = () => ({ left: 0, top: 0, width: 320, height: 240 });
  return canvas;
}

describe('x-desktop-rotation-input', () => {
  test('touch controls do not remove ArcRotate pointer inputs', () => {
    const canvas = makeCanvas();
    const scene = { pointerX: 0, pointerY: 0 };
    const camera = {
      attachControl: jest.fn(),
      detachControl: jest.fn(),
      inputs: {
        removeByType: jest.fn(),
        remove: jest.fn(),
        attached: { pointers: {} },
      },
    };

    const disposer = installTouchControls({ canvas, scene, camera });
    expect(disposer).toBeTruthy();

    expect(camera.inputs.removeByType).not.toHaveBeenCalled();
    expect(camera.inputs.remove).not.toHaveBeenCalled();
  });
});

