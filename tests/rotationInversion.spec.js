import { applyIncrement, identityQ, rotateForward } from '../public/vr/rotation-core.js';

// Simple camera mock: right axis changes when we yaw the accumulated quaternion; here we simulate dynamic right axis by rotating base axes.
function rotateVec(q,v){
  // reuse rotation-core internal math via rotateForward pattern; reimplement minimal
  const vx=v[0], vy=v[1], vz=v[2];
  const {x:qx,y:qy,z:qz,w:qw}=q;
  const tx=2*(qy*vz - qz*vy);
  const ty=2*(qz*vx - qx*vz);
  const tz=2*(qx*vy - qy*vx);
  return [
    vx + qw*tx + (qy*tz - qz*ty),
    vy + qw*ty + (qz*tx - qx*tz),
    vz + qw*tz + (qx*ty - qy*tx)
  ];
}

function quatMul(a,b){
  return { w: a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
           x: a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
           y: a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
           z: a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w };
}
function quatFromAxisAngle(axis, angle){
  const h=angle*0.5; const s=Math.sin(h); return { x:axis[0]*s, y:axis[1]*s, z:axis[2]*s, w:Math.cos(h) };
}

function yawQuat(angle){ return quatFromAxisAngle([0,1,0], angle); }

// Derive camera right from accumulator (rotate world X by accQ)
function cameraRightFrom(accQ){ return rotateVec(accQ, [1,0,0]); }

// Apply a series of yaw rotations to simulate user rotating molecule 180 degrees, then test vertical pitch direction.

describe('rotation inversion regression', () => {
  test('pitch direction consistent after 180 yaw', () => {
    let acc = identityQ();
    // initial camera right (world X)
    let right = cameraRightFrom(acc);
    // Simulate a small positive pitch input (drag up) before yaw inversion
    const before = applyIncrement(acc, 0, 0.02, right, { sens:1, maxStep:0.1 });
    const forwardBefore = rotateForward(before);
    // Expect forwardBefore.y < 0 (since positive pitch up usually makes forward tilt down in right-handed system)
    const signBefore = Math.sign(forwardBefore[1]);

    // Apply yaw ~ PI (rotate molecule)
    const steps = 20; const totalYaw = Math.PI; for(let i=0;i<steps;i++){ acc = applyIncrement(acc, (totalYaw/steps), 0, cameraRightFrom(acc), { sens:1, maxStep:0.2 }); }

    // Now apply same pitch increment after yaw inversion
    right = cameraRightFrom(acc);
    const after = applyIncrement(acc, 0, 0.02, right, { sens:1, maxStep:0.1 });
    const forwardAfter = rotateForward(after);
    const signAfter = Math.sign(forwardAfter[1]);

    // Signs should match (no inversion of vertical response)
    expect(signAfter).toBe(signBefore);
  });
});
