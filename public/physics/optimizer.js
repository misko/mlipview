// Geometry (and future cell) optimization core.
// Provides createOptimizer(state, { forceProvider, settings }) returning step() and runUntil().
// Supports algorithms: 'sd' (steepest descent) and 'cg' (Fletcher-Reeves conjugate gradient).

export function createOptimizer(simState, { forceProvider, settings = {} }) {
  const opts = {
    algorithm: 'sd', // 'sd' | 'cg'
    maxSteps: 500,
    forceTol: 1e-3,
    stressTol: 5e-3,
    gamma: 0.05, // SD step size scaling
    cgReset: 10,
    logEvery: 25,
    wantStress: true,
    ...settings,
  };
  const n = simState.Z.length;
  let prevGradient = null; // for CG (flattened force vector)
  let direction = null; // search direction vector

  function flattenForces(forces) {
    const out = new Float64Array(3 * n);
    for (let i = 0; i < n; i++) {
      const f = forces[i];
      out[3 * i] = f[0];
      out[3 * i + 1] = f[1];
      out[3 * i + 2] = f[2];
    }
    return out;
  }

  function maxForceMag(flatF) {
    let m = 0;
    for (let i = 0; i < flatF.length; i += 3) {
      const fx = flatF[i], fy = flatF[i + 1], fz = flatF[i + 2];
      const mag = Math.hypot(fx, fy, fz);
      if (mag > m) m = mag;
    }
    return m;
  }

  async function step() {
    const pos = simState.pos;
    // Build positions array for provider
    const coords = [];
    for (let i = 0; i < n; i++) coords.push([pos[3 * i], pos[3 * i + 1], pos[3 * i + 2]]);
    const cell = simState.box ? [
      [simState.box.a.x, simState.box.a.y, simState.box.a.z],
      [simState.box.b.x, simState.box.b.y, simState.box.b.z],
      [simState.box.c.x, simState.box.c.y, simState.box.c.z],
    ] : null;
    const { energy, forces, stress } = await forceProvider.compute({ elements: Array.from(simState.Z), positions: coords, cell, wantStress: opts.wantStress });
    // Flatten forces (negative gradient of energy w.r.t positions)
    const flatF = flattenForces(forces);
    const maxF = maxForceMag(flatF);
    let converged = maxF < opts.forceTol;
    // --- Position update ---
    if (opts.algorithm === 'sd') {
      for (let i = 0; i < n; i++) {
        pos[3 * i] += opts.gamma * flatF[3 * i];
        pos[3 * i + 1] += opts.gamma * flatF[3 * i + 1];
        pos[3 * i + 2] += opts.gamma * flatF[3 * i + 2];
      }
      prevGradient = flatF;
    } else if (opts.algorithm === 'cg') {
      if (!prevGradient) {
        prevGradient = flatF;
        direction = flatF.slice();
      } else {
        // Fletcher-Reeves beta = (g_k^T g_k) / (g_{k-1}^T g_{k-1}) where g = -F; here F already acts like -grad, so use F magnitude.
        let num = 0, denom = 0;
        for (let i = 0; i < flatF.length; i++) { num += flatF[i] * flatF[i]; denom += prevGradient[i] * prevGradient[i]; }
        const beta = denom > 0 ? (num / denom) : 0;
        for (let i = 0; i < flatF.length; i++) {
          direction[i] = flatF[i] + beta * direction[i];
        }
        prevGradient = flatF;
        // Periodic reset to avoid drift in noisy landscapes
        if (opts.cgReset && Math.random() < 1 / opts.cgReset) {
          direction = flatF.slice();
        }
      }
      // Simple fixed step length for now; could add line search later.
      for (let i = 0; i < n; i++) {
        pos[3 * i] += opts.gamma * direction[3 * i];
        pos[3 * i + 1] += opts.gamma * direction[3 * i + 1];
        pos[3 * i + 2] += opts.gamma * direction[3 * i + 2];
      }
    } else {
      throw new Error('Unknown algorithm: ' + opts.algorithm);
    }
    simState.energy = energy;
    return { energy, maxForce: maxF, converged };
  }

  async function runUntil(onUpdate) {
    for (let stepIdx = 0; stepIdx < opts.maxSteps; stepIdx++) {
      const res = await step();
      if (onUpdate) onUpdate({ step: stepIdx, ...res });
      if (res.converged) return { ...res, steps: stepIdx + 1 };
    }
    return { converged: false };
  }

  return { step, runUntil, options: opts };
}
