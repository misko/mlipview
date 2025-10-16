// Unit test for joystickDelta mapping function used by VR drag

import { joystickDelta } from '../public/vr/spherical-drag-core.js';

describe('joystickDelta', () => {
  test('deadzone yields zero', () => {
    expect(joystickDelta(0.01, { deadzone: 0.08 })).toBe(0);
  });
  test('forward (negative) increases radius', () => {
    const d = joystickDelta(-1.0, { gain: 0.5, deadzone: 0.05 });
    expect(d).toBeCloseTo(0.5, 5);
  });
  test('backward (positive) decreases radius', () => {
    const d = joystickDelta(0.6, { gain: 0.5, deadzone: 0.05 });
    expect(d).toBeCloseTo(-0.3, 5);
  });
});
