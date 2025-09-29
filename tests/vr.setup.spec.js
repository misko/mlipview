import { createVRSupport } from '../public/vr/setup.js';

describe('VR scaffolding', () => {
  test('exposes expected API', () => {
    const scene = {}; // minimal stub
    const vr = createVRSupport(scene, {});
    expect(typeof vr.init).toBe('function');
    expect(typeof vr.enterVR).toBe('function');
    expect(typeof vr.exitVR).toBe('function');
  });
});
