// Variable-cell optimizer: adjusts cell vectors using stress feedback.
// Simple linear strain update: strain = -cellGamma * stress (Voigt) with capping.
// Only uses diagonal + shear components present; if stress absent, it's a no-op.

import { applyStrain, stressVoigt } from './sim-model.js';

export function createCellOptimizer(simState, { forceProvider, settings = {} }) {
  const opts = {
    maxSteps: 200,
    stressTol: 5e-3,
    cellGamma: 1e-3, // strain per unit stress
    maxStrainStep: 5e-3, // cap any component
    isotropic: false, // if true, average hydrostatic component only
    logEvery: 50,
    ...settings,
  };
  function maxAbsStress(voigt) {
    let m = 0; for (let i=0;i<3;i++){ const v=Math.abs(voigt[i]); if(v>m)m=v; }
    // include shear terms (scaled) if present
    for (let i=3;i<6;i++){ const v=Math.abs(voigt[i]); if(v>m)m=v; }
    return m;
  }
  async function step() {
    // Build positions for provider (cell stress may depend on positions even if we ignore forces here)
    const n = simState.Z.length;
    const coords = []; const pos = simState.pos;
    for (let i=0;i<n;i++) coords.push([pos[3*i], pos[3*i+1], pos[3*i+2]]);
    const cell = simState.box ? [
      [simState.box.a.x, simState.box.a.y, simState.box.a.z],
      [simState.box.b.x, simState.box.b.y, simState.box.b.z],
      [simState.box.c.x, simState.box.c.y, simState.box.c.z],
    ] : null;
    const { stress } = await forceProvider.compute({ elements: Array.from(simState.Z), positions: coords, cell, wantStress: true });
    if (!stress || !stress.voigt) {
      return { changed:false, converged:true, maxStress:0 };
    }
    const v = stress.voigt; // [xx,yy,zz,yz,xz,xy]
    const maxS = maxAbsStress(v);
    const converged = maxS < opts.stressTol;
    if (converged) return { changed:false, converged, maxStress:maxS };
    // Determine strain increment
    let strain = new Float64Array(6);
    if (opts.isotropic) {
      const p = (v[0]+v[1]+v[2])/3; // hydrostatic
      const iso = -opts.cellGamma * p;
      strain[0]=strain[1]=strain[2] = clamp(iso, opts.maxStrainStep);
    } else {
      for (let i=0;i<6;i++) {
        strain[i] = clamp(-opts.cellGamma * v[i], opts.maxStrainStep);
      }
    }
    applyStrain(simState, strain);
    return { changed:true, converged:false, maxStress:maxS };
  }
  async function runUntil(onUpdate){
    for (let i=0;i<opts.maxSteps;i++){
      const res = await step();
      if (onUpdate) onUpdate({ step:i, ...res });
      if (res.converged) return { ...res, steps:i+1 };
    }
    return { converged:false };
  }
  return { step, runUntil, options:opts };
}

function clamp(x, limit){ if (x>limit) return limit; if (x<-limit) return -limit; return x; }
