// Shared simulation data model for relaxation & MD
// Units: length Å, energy eV, force eV/Å, stress eV/Å^3, temperature K.

export function createSimulationState({ Z, positions, cell, masses }={}) {
  const n = Z?.length||0;
  const pos = new Float64Array(3*n);
  if(positions){
    for(let i=0;i<n;i++){ const p=positions[i]; pos[3*i]=p[0]; pos[3*i+1]=p[1]; pos[3*i+2]=p[2]; }
  }
  const vel = new Float64Array(3*n); // zero initial velocities
  const forces = new Float64Array(3*n); // updated in-place
  const mass = new Float64Array(n);
  if(masses){ for(let i=0;i<n;i++) mass[i]=masses[i]; } else { for(let i=0;i<n;i++) mass[i]=defaultMassForZ(Z[i]); }
  const stress = { xx:0, yy:0, zz:0, xy:0, xz:0, yz:0 };
  const box = cell ? normalizeCell(cell) : null;
  return {
    units: baseUnits(),
    Z: Int32Array.from(Z||[]),
    pos, vel, forces, mass,
    step: 0,
    energy: null,
    kinetic: 0,
    temperature: 0,
    stress,
    box,
    meta: { created: Date.now() },
    flags: { variableCell: false }
  };
}

export function baseUnits(){
  return { length:'A', energy:'eV', force:'eV/A', stress:'eV/A^3', temperature:'K' };
}

export function defaultMassForZ(Z){
  // Very rough approximate atomic masses (amu) for light elements, fallback 12.
  const table={1:1.008,6:12.011,7:14.007,8:15.999,16:32.06};
  return table[Z]||12.0;
}

export function zeroForces(state){ state.forces.fill(0); }

export function stressNorm(s){ return Math.sqrt(s.xx*s.xx + s.yy*s.yy + s.zz*s.zz + 2*(s.xy*s.xy + s.xz*s.xz + s.yz*s.yz)); }

export function maxForce(state){
  let m=0; const f=state.forces; for(let i=0;i<f.length;i+=3){ const fx=f[i],fy=f[i+1],fz=f[i+2]; const mag=Math.hypot(fx,fy,fz); if(mag>m) m=mag; }
  return m;
}

export function rmsForce(state){
  const f=state.forces; let sum=0; let n=0; for(let i=0;i<f.length;i+=3){ const fx=f[i],fy=f[i+1],fz=f[i+2]; sum+=fx*fx+fy*fy+fz*fz; n++; }
  return n? Math.sqrt(sum/n) : 0;
}

export function computeKineticAndTemp(state){
  const kB_eV_per_K = 8.617333262145e-5; // Boltzmann constant in eV/K
  const vel=state.vel, mass=state.mass; let K=0; for(let i=0;i<mass.length;i++){ const m=mass[i]; const vx=vel[3*i], vy=vel[3*i+1], vz=vel[3*i+2]; // mass in amu -> convert to kg? We keep internal K relative scale.
    // For now treat masses in arbitrary amu and compute pseudo kinetic energy in eV scaling constant factor; refine later.
    K += 0.5 * m * (vx*vx + vy*vy + vz*vz);
  }
  state.kinetic = K; // pseudo
  const dof = 3*mass.length;
  state.temperature = dof>0 ? (2*K/(dof*kB_eV_per_K)) : 0;
}

export function normalizeCell(cell){
  // Accept 3x3 array [[ax,ay,az],...]
  if(!Array.isArray(cell) || cell.length!==3) return null;
  const a=cell[0], b=cell[1], c=cell[2];
  const mat = [ ...a, ...b, ...c ];
  const vol = cellVolume(a,b,c);
  const lengths=[vecLen(a),vecLen(b),vecLen(c)];
  const angles=[ angleDeg(b,c), angleDeg(a,c), angleDeg(a,b) ]; // alpha=bc, beta=ac, gamma=ab
  return { a:{x:a[0],y:a[1],z:a[2]}, b:{x:b[0],y:b[1],z:b[2]}, c:{x:c[0],y:c[1],z:c[2]}, matrix:mat, volume:vol, lengths, angles };
}

function vecLen(v){ return Math.hypot(v[0],v[1],v[2]); }
function dot(a,b){ return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }
function angleDeg(a,b){ const d=dot(a,b)/(vecLen(a)*vecLen(b)||1); return Math.acos(Math.max(-1,Math.min(1,d)))*180/Math.PI; }
function cellVolume(a,b,c){
  // triple product |a · (b × c)|
  const bx=b[0],by=b[1],bz=b[2];
  const cx=c[0],cy=c[1],cz=c[2];
  const cross=[ by*cz - bz*cy, bz*cx - bx*cz, bx*cy - by*cx ];
  return Math.abs(a[0]*cross[0] + a[1]*cross[1] + a[2]*cross[2]);
}

export function applyStrain(state, strain){
  // strain: symmetric 3x3 in Voigt order [e_xx,e_yy,e_zz,e_yz,e_xz,e_xy]
  if(!state.box) return;
  const { a,b,c } = state.box; const A=[a.x,a.y,a.z], B=[b.x,b.y,b.z], C=[c.x,c.y,c.z];
  const S = strain; // length 6
  function deform(v){
    const x=v[0],y=v[1],z=v[2];
    return [
      x + S[0]*x + S[5]*y + S[4]*z,
      y + S[5]*x + S[1]*y + S[3]*z,
      z + S[4]*x + S[3]*y + S[2]*z,
    ];
  }
  const A2=deform(A), B2=deform(B), C2=deform(C);
  state.box = normalizeCell([A2,B2,C2]);
}

export function stressVoigt(stress){ return [stress.xx, stress.yy, stress.zz, stress.yz, stress.xz, stress.xy]; }
export function setStressFromVoigt(stress, arr){ [stress.xx,stress.yy,stress.zz,stress.yz,stress.xz,stress.xy] = arr; }

export function updatePositions(state, newPos){ for(let i=0;i<state.Z.length;i++){ state.pos[3*i]=newPos[3*i]; state.pos[3*i+1]=newPos[3*i+1]; state.pos[3*i+2]=newPos[3*i+2]; } }

// Placeholder hook for unit conversions (if FairChem differs later)
export function ensureUnits(state, fairchemMeta){ return state; }
