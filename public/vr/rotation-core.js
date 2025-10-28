// Pure rotation core (framework agnostic)
// Quaternions represented as {x,y,z,w}

function quatMul(a, b) {
  return {
    w: a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z,
    x: a.w * b.x + a.x * b.w + a.y * b.z - a.z * b.y,
    y: a.w * b.y - a.x * b.z + a.y * b.w + a.z * b.x,
    z: a.w * b.z + a.x * b.y - a.y * b.x + a.z * b.w,
  };
}
function quatNormalize(q) {
  const len = Math.hypot(q.w, q.x, q.y, q.z) || 1;
  return { x: q.x / len, y: q.y / len, z: q.z / len, w: q.w / len };
}
function quatFromAxisAngle(axis, angle) {
  const h = angle * 0.5;
  const s = Math.sin(h);
  return { x: axis[0] * s, y: axis[1] * s, z: axis[2] * s, w: Math.cos(h) };
}
function rotateVec(q, v) {
  // v:[x,y,z]
  // q * v * q^-1 using quaternion math (optimized)
  const vx = v[0],
    vy = v[1],
    vz = v[2];
  const qx = q.x,
    qy = q.y,
    qz = q.z,
    qw = q.w;
  // t = 2 * cross(q.xyz, v)
  const tx = 2 * (qy * vz - qz * vy);
  const ty = 2 * (qz * vx - qx * vz);
  const tz = 2 * (qx * vy - qy * vx);
  // v' = v + qw*t + cross(q.xyz, t)
  return [
    vx + qw * tx + (qy * tz - qz * ty),
    vy + qw * ty + (qz * tx - qx * tz),
    vz + qw * tz + (qx * ty - qy * tx),
  ];
}

// Apply incremental yaw (world Y) and pitch (camera right axis) to accumulated quaternion.
// dYaw/dPitch are raw small deltas (already sign-adjusted by caller). Sensitivity & clamps applied here.
import { __count } from '../util/funcCount.js';
export function applyIncrement(accQ, dYaw, dPitch, cameraRight, opts = {}) {
  __count('vrRotation#applyIncrement');
  const sens = opts.sens != null ? opts.sens : 0.9;
  const maxStep = opts.maxStep != null ? opts.maxStep : 0.05; // radians per frame
  const clamp = (v) => Math.max(-maxStep, Math.min(maxStep, v));
  const cy = clamp(dYaw) * sens;
  const cp = clamp(dPitch) * sens;
  if (Math.abs(cy) < 1e-8 && Math.abs(cp) < 1e-8) return accQ;
  const qYaw = quatFromAxisAngle([0, 1, 0], cy);
  const rightNormLen = Math.hypot(cameraRight[0], cameraRight[1], cameraRight[2]) || 1;
  const rAxis = [
    cameraRight[0] / rightNormLen,
    cameraRight[1] / rightNormLen,
    cameraRight[2] / rightNormLen,
  ];
  const qPitch = quatFromAxisAngle(rAxis, cp);
  // Order: apply yaw then pitch in world/user frame -> dq = qYaw * qPitch
  const dq = quatMul(qYaw, qPitch);
  const out = quatNormalize(quatMul(dq, accQ));
  return out;
}

export function identityQ() {
  __count('vrRotation#identityQ');
  return { x: 0, y: 0, z: 0, w: 1 };
}
export function rotateForward(q) {
  __count('vrRotation#rotateForward');
  return rotateVec(q, [0, 0, 1]);
}

// For tests
export const __private = { quatMul, quatFromAxisAngle, rotateVec, quatNormalize };
