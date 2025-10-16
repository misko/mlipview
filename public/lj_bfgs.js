// Minimal Lennard-Jones potential and BFGS optimizer matching ASE defaults.
// Parity target: initial water energy 2.802523 -> relaxed energy -2.983548 (fmax 0.05)
// Assumptions: epsilon=1.0, sigma=1.0 global; rc = 3*sigma with energy shift (no smooth cutoff)

// LJ with standard energy shift at cutoff so U(rc)=0.
// Parity target: water initial energy 2.802523 -> relaxed -2.983548 using rc=3*sigma, epsilon=1, sigma=1.
export function ljEnergyForces(positions, { epsilon = 1.0, sigma = 1.0 } = {}) {
  const n = positions.length; // positions: [ [x,y,z], ... ]
  const forces = Array.from({ length: n }, () => [0, 0, 0]);
  let energy = 0;
  const rc = 3 * sigma;
  const invrc = 1 / rc;
  const src = sigma * invrc; // sigma/rc
  const src2 = src * src;
  const src6 = src2 * src2 * src2;
  const src12 = src6 * src6;
  const e0 = 4 * epsilon * (src12 - src6); // U(rc) (negative); shift so potential zero at cutoff
  for (let i = 0; i < n; i++) {
    const pi = positions[i];
    for (let j = i + 1; j < n; j++) {
      const pj = positions[j];
      const dx = pj[0] - pi[0];
      const dy = pj[1] - pi[1];
      const dz = pj[2] - pi[2];
      const r2 = dx * dx + dy * dy + dz * dz;
      if (r2 === 0) continue;
      const r = Math.sqrt(r2);
      if (r > rc) continue; // truncated
      const invr = 1 / r;
      const sr = sigma * invr;
      const sr2 = sr * sr;
      const sr6 = sr2 * sr2 * sr2;
      const sr12 = sr6 * sr6;
  let pairE = 4 * epsilon * (sr12 - sr6) - e0; // shifted so U(rc)=0 (adds ~0.00548 per pair at rc=3)
      energy += pairE;
      // force magnitude factor: dU/dr times direction
      // dU/dr = 24*epsilon*(2*sr12 - sr6)/r
  // Standard pair force magnitude along vector from i->j is F = 24*epsilon*(2*sr12 - sr6)/r
  const forceMagOverR = 24 * epsilon * (2 * sr12 - sr6) * invr * invr; // (1/r)*dU/dr but sign already accounts
  // Force on atom i due to j is -forceMag * (r_j - r_i)/r -> negative of below
  const fx = forceMagOverR * dx;
  const fy = forceMagOverR * dy;
  const fz = forceMagOverR * dz;
  // Apply opposite signs: i gets -f, j gets +f
  forces[i][0] -= fx; forces[i][1] -= fy; forces[i][2] -= fz;
  forces[j][0] += fx; forces[j][1] += fy; forces[j][2] += fz;
    }
  }
  return { energy, forces };
}

// Basic BFGS optimizer with line search (Armijo backtracking).
export async function bfgsOptimize({ positions, fmax = 0.05, maxSteps = 1000, maxStep = 0.2, compute }) {
  // positions is mutable array of [x,y,z]
  const n = positions.length;
  const dim = 3 * n;
  function flatten(posArr) {
    const out = new Float64Array(dim);
    for (let i = 0; i < n; i++) { const p = posArr[i]; out[3 * i] = p[0]; out[3 * i + 1] = p[1]; out[3 * i + 2] = p[2]; }
    return out;
  }
  function unflatten(vec) {
    for (let i = 0; i < n; i++) { positions[i][0] = vec[3 * i]; positions[i][1] = vec[3 * i + 1]; positions[i][2] = vec[3 * i + 2]; }
  }
  function flattenForces(forces) {
    const g = new Float64Array(dim);
    for (let i = 0; i < n; i++) { const f = forces[i]; g[3 * i] = -f[0]; g[3 * i + 1] = -f[1]; g[3 * i + 2] = -f[2]; } // gradient = -forces
    return g;
  }
  function maxForceMag(forces) {
    let m = 0; for (const f of forces) { const mag = Math.hypot(f[0], f[1], f[2]); if (mag > m) m = mag; } return m;
  }
  // Initial evaluation
  let { energy: E, forces: F } = await compute(positions);
  let g = flattenForces(F);
  let x = flatten(positions);
  let H = identity(dim); // inverse Hessian approximation
  const history = [{ step: 0, energy: E, fmax: maxForceMag(F) }];
  if (history[0].fmax < fmax) return { converged: true, energy: E, steps: 0, history };
  for (let k = 0; k < maxSteps; k++) {
    // p = -H g
    const p = new Float64Array(dim);
    for (let i = 0; i < dim; i++) { let s = 0; const row = H[i]; for (let j = 0; j < dim; j++) s += row[j] * g[j]; p[i] = -s; }
    // Limit maximum component displacement
    let maxComp = 0; for (let i = 0; i < dim; i++) { const a = Math.abs(p[i]); if (a > maxComp) maxComp = a; }
    if (maxComp > maxStep) {
      const scale = maxStep / maxComp; for (let i = 0; i < dim; i++) p[i] *= scale;
    }
    // Line search (Armijo)
    const c1 = 1e-4; let alpha = 1.0; const gdotp = dot(g, p);
    const x_old = x.slice(); let newE, newF, newG;
    while (true) {
      for (let i = 0; i < dim; i++) x[i] = x_old[i] + alpha * p[i];
      unflatten(x);
      ({ energy: newE, forces: newF } = await compute(positions));
      if (newE <= E + c1 * alpha * gdotp || alpha < 1e-6) {
        newG = flattenForces(newF);
        break;
      }
      alpha *= 0.5;
    }
    // BFGS update
    const s = new Float64Array(dim); for (let i = 0; i < dim; i++) s[i] = x[i] - x_old[i];
    const y = new Float64Array(dim); for (let i = 0; i < dim; i++) y[i] = newG[i] - g[i];
    const ys = dot(y, s);
    if (ys > 1e-12) {
      const Hy = matVec(H, y);
      const yHy = dot(y, Hy);
      // H_{k+1} = (I - s y^T / ys) H (I - y s^T / ys) + (s s^T)/ys
      // Implement via outer products (dense, small dim ok for tiny systems)
      const rho = 1 / ys;
      // Compute (I - rho s y^T) H (I - rho y s^T)
      const temp = Array.from({ length: dim }, () => new Float64Array(dim));
      for (let i = 0; i < dim; i++) {
        for (let j = 0; j < dim; j++) {
          // First product H (I - rho y s^T)
          let Hij = H[i][j];
          let corr = 0;
          for (let q = 0; q < dim; q++) corr += H[i][q] * (-rho * y[q] * s[j]);
          temp[i][j] = Hij + corr;
        }
      }
      const newH = Array.from({ length: dim }, () => new Float64Array(dim));
      for (let i = 0; i < dim; i++) {
        for (let j = 0; j < dim; j++) {
          // (I - rho s y^T) temp
            let val = temp[i][j];
            let corr = 0;
            for (let q = 0; q < dim; q++) corr += (-rho * s[i] * y[q]) * temp[q][j];
            newH[i][j] = val + corr + rho * s[i] * s[j];
        }
      }
      H = newH;
    } else {
      // Reset if update not positive definite
      H = identity(dim);
    }
    E = newE; F = newF; g = newG;
    const fmaxNow = maxForceMag(F);
    history.push({ step: k + 1, energy: E, fmax: fmaxNow });
    if (fmaxNow < fmax) return { converged: true, energy: E, steps: k + 1, history };
  }
  return { converged: false, energy: E, steps: maxSteps, history };
}

function identity(n) { return Array.from({ length: n }, (_, i) => { const row = new Float64Array(n); row[i] = 1; return row; }); }
function dot(a, b) { let s = 0; for (let i = 0; i < a.length; i++) s += a[i] * b[i]; return s; }
function matVec(M, v) { const out = new Float64Array(v.length); for (let i = 0; i < M.length; i++) { let s = 0; const row = M[i]; for (let j = 0; j < v.length; j++) s += row[j] * v[j]; out[i] = s; } return out; }

// Convenience: run water optimization
export async function optimizeWaterExample() {
  const d = 0.9575;
  const t = Math.PI / 180 * 104.51;
  const positions = [
    [d, 0, 0],
    [d * Math.cos(t), d * Math.sin(t), 0],
    [0, 0, 0]
  ];
  const initial = ljEnergyForces(positions);
  const res = await bfgsOptimize({ positions, fmax: 0.05, maxSteps: 500, maxStep: 0.2, compute: async (p) => ljEnergyForces(p) });
  return { initial, final: res, positions };
}

// Incremental BFGS stepper: returns object with step() performing one iteration and exposing trace.
export function createBFGSStepper({ positions, fmax = 0.05, maxStep = 0.2, compute }) {
  const n = positions.length; const dim = 3 * n;
  const history = [];
  let H = identity(dim);
  let x = flatten(positions);
  let evalCache = null;
  async function ensureEval(){ if(!evalCache){ const { energy, forces } = await compute(positions); evalCache = { energy, forces, g: flattenForces(forces), fmax: maxForceMag(forces) }; history.push({ step: history.length, energy, fmax: evalCache.fmax }); } return evalCache; }
  function flatten(posArr){ const out=new Float64Array(dim); for(let i=0;i<n;i++){ const p=posArr[i]; out[3*i]=p[0]; out[3*i+1]=p[1]; out[3*i+2]=p[2]; } return out; }
  function unflatten(vec){ for(let i=0;i<n;i++){ positions[i][0]=vec[3*i]; positions[i][1]=vec[3*i+1]; positions[i][2]=vec[3*i+2]; } }
  function flattenForces(forces){ const g=new Float64Array(dim); for(let i=0;i<n;i++){ const f=forces[i]; g[3*i]=-f[0]; g[3*i+1]=-f[1]; g[3*i+2]=-f[2]; } return g; }
  function maxForceMag(forces){ let m=0; for(const f of forces){ const mag=Math.hypot(f[0],f[1],f[2]); if(mag>m)m=mag; } return m; }
  function dot(a,b){ let s=0; for(let i=0;i<a.length;i++) s+=a[i]*b[i]; return s; }
  return {
    history,
    // Direct reference to mutable positions array (array of [x,y,z]) so caller can sync to state
    positions,
    async step(){
      const { energy:E, forces:F, g } = await ensureEval();
      if (history.length>0 && history[history.length-1].energy!==E) {
        // already recorded
      }
      if (Math.max(...F.map(f=>Math.hypot(...f))) < fmax) return { converged:true, energy:E, fmax: history[history.length-1].fmax };
      // p = -H g
      const dimLocal = g.length; const p = new Float64Array(dimLocal);
      for (let i=0;i<dimLocal;i++){ const row=H[i]; let s=0; for(let j=0;j<dimLocal;j++) s+=row[j]*g[j]; p[i]=-s; }
      // limit step
      let maxComp=0; for(let i=0;i<p.length;i++){ const a=Math.abs(p[i]); if(a>maxComp) maxComp=a; }
      if(maxComp>maxStep){ const scale=maxStep/maxComp; for(let i=0;i<p.length;i++) p[i]*=scale; }
      // line search
      const c1=1e-4; let alpha=1.0; const gdotp = dot(g,p); const x_old = x.slice(); let newEval=null;
      while(true){
        for(let i=0;i<x.length;i++) x[i]=x_old[i] + alpha*p[i];
        unflatten(x);
        const { energy, forces } = await compute(positions);
        const newG = flattenForces(forces);
        if (energy <= E + c1*alpha*gdotp || alpha < 1e-6) { newEval = { energy, forces, g:newG }; break; }
        alpha *= 0.5;
      }
      // BFGS update
      const s = new Float64Array(p.length); for(let i=0;i<s.length;i++) s[i]=alpha*p[i];
      const y = new Float64Array(p.length); for(let i=0;i<y.length;i++) y[i]=newEval.g[i]-g[i];
      const ys = dot(y,s);
      if (ys>1e-12){
        const rho = 1/ys;
        // (I - rho s y^T) H (I - rho y s^T) + rho s s^T
        const dimL = p.length; const temp = Array.from({length:dimL},()=>new Float64Array(dimL));
        for(let i=0;i<dimL;i++){
          for(let j=0;j<dimL;j++){
            let Hij=H[i][j]; let corr=0; for(let q=0;q<dimL;q++) corr += H[i][q]*(-rho*y[q]*s[j]);
            temp[i][j]=Hij + corr;
          }
        }
        const newH = Array.from({length:dimL},()=>new Float64Array(dimL));
        for(let i=0;i<dimL;i++){
          for(let j=0;j<dimL;j++){
            let val=temp[i][j]; let corr=0; for(let q=0;q<dimL;q++) corr += (-rho*s[i]*y[q])*temp[q][j];
            newH[i][j]= val + corr + rho*s[i]*s[j];
          }
        }
        H = newH;
      } else {
        H = identity(p.length);
      }
      evalCache = { ...newEval, fmax: maxForceMag(newEval.forces) };
      history.push({ step: history.length, energy: evalCache.energy, fmax: evalCache.fmax });
      return { converged:false, energy: evalCache.energy, fmax: evalCache.fmax };
    },
    getCurrent: async ()=>{ const e=await ensureEval(); return { energy:e.energy, forces:e.forces, fmax:e.fmax }; }
  };
}
