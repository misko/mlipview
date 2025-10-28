import { applyIncrement, identityQ, rotateForward } from '../public/vr/rotation-core.js';

function rotateVec(q, v) {
  const vx = v[0],
    vy = v[1],
    vz = v[2];
  const { x: qx, y: qy, z: qz, w: qw } = q;
  const tx = 2 * (qy * vz - qz * vy);
  const ty = 2 * (qz * vx - qx * vz);
  const tz = 2 * (qx * vy - qy * vx);
  return [
    vx + qw * tx + (qy * tz - qz * ty),
    vy + qw * ty + (qz * tx - qx * tz),
    vz + qw * tz + (qx * ty - qy * tx),
  ];
}

function cameraRightFrom(accQ) {
  return rotateVec(accQ, [1, 0, 0]);
}

describe('x-rotation-inversion', () => {
  test('pitch direction remains stable after 180° yaw', () => {
    let acc = identityQ();
    const rightInitial = cameraRightFrom(acc);
    const beforePitch = applyIncrement(acc, 0, 0.02, rightInitial, { sens: 1, maxStep: 0.1 });
    const forwardBefore = rotateForward(beforePitch);
    const signBefore = Math.sign(forwardBefore[1]);

    const steps = 24;
    const totalYaw = Math.PI;
    for (let i = 0; i < steps; i++) {
      const stepYaw = totalYaw / steps;
      const rightAxis = cameraRightFrom(acc);
      acc = applyIncrement(acc, stepYaw, 0, rightAxis, { sens: 1, maxStep: 0.2 });
    }

    const rightAfter = cameraRightFrom(acc);
    const afterPitch = applyIncrement(acc, 0, 0.02, rightAfter, { sens: 1, maxStep: 0.1 });
    const forwardAfter = rotateForward(afterPitch);
    const signAfter = Math.sign(forwardAfter[1]);

    expect(signAfter).toBe(signBefore);
  });
});
