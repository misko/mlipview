export function createForceField(molState, { kBond = 200, r0 = 1.1, ljEpsilon = 0.02, ljSigma = 3.4 } = {}) {
  function approximateVolume() {
    const c = molState.cell;
    if (c && c.enabled && c.a && c.b && c.c) {
      const ax=c.a.x, ay=c.a.y, az=c.a.z;
      const bx=c.b.x, by=c.b.y, bz=c.b.z;
      const cx=c.c.x, cy=c.c.y, cz=c.c.z;
      const vx = ay*bz - az*by;
      const vy = az*bx - ax*bz;
      const vz = ax*by - ay*bx;
      return Math.abs(vx*cx + vy*cy + vz*cz) || 1;
    }
    let minX=Infinity,minY=Infinity,minZ=Infinity,maxX=-Infinity,maxY=-Infinity,maxZ=-Infinity;
    for (const p of molState.positions) { if(p.x<minX)minX=p.x; if(p.y<minY)minY=p.y; if(p.z<minZ)minZ=p.z; if(p.x>maxX)maxX=p.x; if(p.y>maxY)maxY=p.y; if(p.z>maxZ)maxZ=p.z; }
    return (maxX-minX+1e-6)*(maxY-minY+1e-6)*(maxZ-minZ+1e-6);
  }
  function computeForces() {
    const N = molState.positions.length;
    let F = molState.dynamics.forces;
    let energy = 0;
    // Ensure forces array is initialized and zeroed
    if (!Array.isArray(F) || F.length !== N) {
      F = molState.dynamics.forces = Array.from({ length:N }, ()=>({ x:0,y:0,z:0 }));
    } else {
      for (let i=0;i<N;i++){ F[i].x=0; F[i].y=0; F[i].z=0; }
    }
    const forces = F;
    if (!molState.dynamics.stress) {
      molState.dynamics.stress = { xx:0,yy:0,zz:0,xy:0,xz:0,yz:0 };
    }
    for (const b of molState.bonds) {
      const a = molState.positions[b.i];
      const c = molState.positions[b.j];
      const dx = c.x-a.x, dy=c.y-a.y, dz=c.z-a.z;
      const r = Math.sqrt(dx*dx+dy*dy+dz*dz) || 1e-9;
      const dr = r - r0;
      const fmag = -kBond * dr;
      const fx = fmag * dx / r, fy = fmag * dy / r, fz = fmag * dz / r;
      forces[b.i].x += fx; forces[b.i].y += fy; forces[b.i].z += fz;
      forces[b.j].x -= fx; forces[b.j].y -= fy; forces[b.j].z -= fz;
      energy += 0.5 * kBond * dr * dr;
    }
    for (let i=0;i<N;i++) {
      const pi = molState.positions[i];
      for (let j=i+1;j<N;j++) {
        const pb = molState.positions[j];
        const dx = pb.x-pi.x, dy=pb.y-pi.y, dz=pb.z-pi.z;
        const r2 = dx*dx+dy*dy+dz*dz;
        if (r2 < 1e-6) continue;
        const invR2 = 1/r2; const sig2 = ljSigma*ljSigma; const sr2 = sig2*invR2;
        const sr6 = sr2*sr2*sr2; const sr12 = sr6*sr6;
        const lj = 4*ljEpsilon*(sr12 - sr6);
        energy += lj;
        const fmag = 24*ljEpsilon*(2*sr12 - sr6) * Math.sqrt(invR2);
        const fx = fmag*dx, fy=fmag*dy, fz=fmag*dz;
        forces[i].x -= fx; forces[i].y -= fy; forces[i].z -= fz;
        forces[j].x += fx; forces[j].y += fy; forces[j].z += fz;
      }
    }
    let Wxx=0,Wyy=0,Wzz=0,Wxy=0,Wxz=0,Wyz=0;
    for (let i=0;i<N;i++) {
      const r = molState.positions[i];
      const f = forces[i];
      Wxx += r.x*f.x; Wyy += r.y*f.y; Wzz += r.z*f.z;
      Wxy += r.x*f.y; Wxz += r.x*f.z; Wyz += r.y*f.z;
    }
    const V = approximateVolume();
    const sig = molState.dynamics.stress;
    sig.xx = -(Wxx)/V; sig.yy = -(Wyy)/V; sig.zz = -(Wzz)/V;
    sig.xy = -(0.5*(Wxy + (Wxy)));
    sig.xz = -(0.5*(Wxz + (Wxz)));
    sig.yz = -(0.5*(Wyz + (Wyz)));
    // Persist potential energy
    if (!molState.dynamics) molState.dynamics = {};
    molState.dynamics.energy = energy;
    return energy;
  }
  return { computeForces };
}
