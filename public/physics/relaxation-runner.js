// Relaxation runner: orchestrates atomic + optional cell optimization.
// Uses existing atomic optimizer (optimizer.js) and cell optimizer (cell-optimizer.js).

import { createOptimizer } from './optimizer.js';
import { createCellOptimizer } from './cell-optimizer.js';

export function createRelaxationRunner(simState, { forceProvider, settings = {} }) {
  const opts = {
    atomicAlgorithm: 'cg',
    atomicEvery: 1, // atomic step frequency
  cellEvery: 2, // attempt cell step every N atomic steps (faster cell convergence)
    maxOuter: 500,
    forceTol: 1e-3,
    stressTol: 5e-3,
    gamma: 0.2,
  cellGamma: 5e-2, // stronger cell update (mock friendly)
  maxStrainStep: 0.2,
    onUpdate: null,
    isotropicCell: true,
    ...settings,
  };
  const atomicOpt = createOptimizer(simState, { forceProvider, settings:{ algorithm: opts.atomicAlgorithm, forceTol: opts.forceTol, gamma: opts.gamma, maxSteps:1 } });
  const cellOpt = createCellOptimizer(simState, { forceProvider, settings:{ stressTol: opts.stressTol, cellGamma: opts.cellGamma, isotropic: opts.isotropicCell, maxSteps:1, maxStrainStep: opts.maxStrainStep } });
  let aborted = false;
  function abort(){ aborted = true; }
  async function run(){
    aborted = false;
    for (let outer=0; outer<opts.maxOuter; outer++) {
      if (aborted) return { aborted:true };
      // Atomic step(s)
      const aRes = await atomicOpt.step();
      let cRes = null;
      if (outer % opts.cellEvery === 0) {
        cRes = await cellOpt.step();
      }
      const forceConv = aRes.converged;
      const stressConv = cRes ? cRes.converged : true; // if no cell step, treat as converged for stress
      if (opts.onUpdate) opts.onUpdate({ step:outer, forceConv, stressConv, maxForce:aRes.maxForce, energy:aRes.energy, maxStress: cRes ? cRes.maxStress : null });
      // Require both force & stress convergence, but allow force to reconverge if cell keeps changing
      if (forceConv && stressConv) return { converged:true, steps:outer+1, energy:aRes.energy };
    }
    return { converged:false };
  }
  return { run, abort, options:opts };
}
