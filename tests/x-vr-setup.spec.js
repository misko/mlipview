import { createVRSupport } from '../public/vr/setup.js';

describe('x-vr setup scaffolding', () => {
  test('createVRSupport exposes init/enter/exit', () => {
    const scene = {
      onBeforeRenderObservable: { add() {} },
      getEngine: () => ({ getRenderingCanvas: () => null }),
    };
    const vr = createVRSupport(scene, {});
    expect(typeof vr.init).toBe('function');
    expect(typeof vr.enterVR).toBe('function');
    expect(typeof vr.exitVR).toBe('function');
  });
});
