import { joystickDelta } from '../public/vr/spherical-drag-core.js';

describe('joystickDelta', () => {
  test('applies deadzone and scaling', () => {
    expect(joystickDelta(0.04)).toBe(0);
    const forward = joystickDelta(-0.8, { gain: 0.5, deadzone: 0.2 });
    expect(forward).toBeGreaterThan(0);
    const backward = joystickDelta(0.6, { gain: 0.5, deadzone: 0.2 });
    expect(backward).toBeLessThan(0);
  });
});
