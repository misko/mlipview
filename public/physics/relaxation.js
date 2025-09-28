// public/physics/relaxation.js
// Real-time molecular geometry relaxation engine

/**
 * Create a real-time relaxation engine for molecular geometry optimization
 * Shows optimization steps live as they happen
 * 
 * @param {Object} options - Configuration options
 * @param {Function} options.getMlip - Function to get current MLIP instance  
 * @param {Object} options.molecule - Molecule object with atoms array
 * @param {Function} options.updateCallback - Called when coordinates are updated
 * @param {Function} options.energyCallback - Called with energy value for plotting
 * @param {Function} options.completionCallback - Called when relaxation completes with status
 * @param {String} options.optimizer - Optimizer type: 'steepest_descent', 'conjugate_gradient'
 * @param {Number} options.maxSteps - Maximum number of optimization steps
 * @param {Number} options.stepSize - Initial step size
 * @param {Number} options.energyTolerance - Energy convergence tolerance
 * @param {Number} options.forceTolerance - Force convergence tolerance (max force component)
 * @param {Number} options.stepDelay - Delay between steps in milliseconds for visualization
 */
export function createRelaxationEngine({
  getMlip,
  molecule,
  updateCallback = () => {},
  energyCallback = () => {},
  optimizer = 'steepest_descent',
  maxSteps = 1000,
  stepSize = 0.01,
  energyTolerance = 1e-4,
  forceTolerance = 0.01,
  stepDelay = 0
} = {}) {

  let isRunning = false;
  let currentStep = 0;
  let lastEnergy = null;
  let lastForces = null;
  let searchDirection = null; // for conjugate gradient
  let beta = 0; // CG parameter
  let stopRequested = false;
  
  const state = {
    step: 0,
    energy: null,
    maxForce: null,
    energyChange: null,
    converged: false,
    status: 'idle' // 'idle', 'running', 'converged', 'max_steps', 'stopped', 'error'
  };

  function reset() {
    isRunning = false;
    stopRequested = false;
    currentStep = 0;
    lastEnergy = null;
    lastForces = null;
    searchDirection = null;
    beta = 0;
    state.step = 0;
    state.energy = null;
    state.maxForce = null;
    state.energyChange = null;
    state.converged = false;
    state.status = 'idle';
  }

  function getAtomPositions() {
    return molecule.atoms.map(atom => ({
      x: atom.pos.x,
      y: atom.pos.y,
      z: atom.pos.z
    }));
  }

  function setAtomPositions(positions) {
    for (let i = 0; i < molecule.atoms.length; i++) {
      molecule.atoms[i].pos.x = positions[i].x;
      molecule.atoms[i].pos.y = positions[i].y;
      molecule.atoms[i].pos.z = positions[i].z;
    }
    // Increment change counter to trigger physics cache invalidation
    if (molecule.changeCounter !== undefined) {
      molecule.changeCounter++;
    }
    updateCallback();
  }

  function vectorLength(vec) {
    return Math.sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
  }

  function vectorAdd(a, b) {
    return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
  }

  function vectorSubtract(a, b) {
    return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
  }

  function vectorScale(vec, scale) {
    return { x: vec.x * scale, y: vec.y * scale, z: vec.z * scale };
  }

  function vectorDot(a, b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }

  function getMaxForce(forces) {
    let maxF = 0;
    for (const force of forces) {
      const fmag = vectorLength(force);
      if (fmag > maxF) maxF = fmag;
    }
    return maxF;
  }

  function babylonForcesToVectors(babylonForces) {
    return babylonForces.map(f => ({ x: f.x, y: f.y, z: f.z }));
  }

  async function computeEnergyAndForces() {
    const mlip = getMlip();
    if (!mlip) throw new Error('No MLIP available');
    // Wrap atoms into primary cell before computing forces (if cell present)
    try {
      if (molecule && molecule.wrapIntoCell) molecule.wrapIntoCell({ refresh: false, bondRefresh: false, log: false });
    } catch {}
    const result = await mlip.compute();
    if (!result || typeof result.energy !== 'number' || !Array.isArray(result.forces)) {
      throw new Error('Invalid MLIP result');
    }
    // Capture stress if present
    let stress = null;
    if (result.stress) {
      stress = result.stress; // could be array length 6, or 3x3; we treat as Voigt length 6 if length==6
    } else if (mlip.lastStress) {
      stress = mlip.lastStress;
    }
    if (stress) {
      try {
        const sxx = stress[0], syy = stress[1], szz = stress[2];
        console.log('[Relaxation][Stress] received from forcefield', {
          hasStress: true,
          length: Array.isArray(stress)? stress.length : 'n/a',
          sxx, syy, szz,
          hydrostatic_p: (Number.isFinite(sxx)&&Number.isFinite(syy)&&Number.isFinite(szz))? Number(((sxx+syy+szz)/3).toFixed(6)) : null
        });
      } catch {}
    } else {
      console.log('[Relaxation][Stress] not present this step');
    }
    
    return {
      energy: result.energy,
      forces: babylonForcesToVectors(result.forces),
      stress
    };
  }

  function tryCellRelax(stress) {
    try {
      if (!stress) return;
      if (!molecule || !molecule.__cellState || !molecule.__cellState.visible) return;
      const cs = molecule.__cellState;
      const { a, b, c } = cs.vectors || {};
      if (!a || !b || !c) return;
      // Interpret stress as Voigt (xx, yy, zz, yz, xz, xy) in some units; we apply only hydrostatic part.
      const sxx = stress[0], syy = stress[1], szz = stress[2];
      if (![sxx,syy,szz].every(v=> Number.isFinite(v))) return;
      const p = (sxx + syy + szz)/3; // simple average
      // If near zero, skip
      if (Math.abs(p) < 1e-4) {
        console.log('[Relaxation][Stress] received (no cell change, |p| small)', { sxx, syy, szz, p: p.toFixed(6) });
        return;
      }
      // Small scaling proportional to pressure; clamp magnitude (heuristic)
      const scaleStep = Math.max(-0.01, Math.min(0.01, -0.002 * p)); // negative p expands, positive compress
      // Cell volume before
      const oldVol = a.x*(b.y*c.z - b.z*c.y) - a.y*(b.x*c.z - b.z*c.x) + a.z*(b.x*c.y - b.y*c.x);
      const newA = a.scale(1+scaleStep);
      const newB = b.scale(1+scaleStep);
      const newC = c.scale(1+scaleStep);
      // Preserve fractional coordinates
      // Build inverse cell matrix with columns a,b,c
      const M = [a.x, b.x, c.x, a.y, b.y, c.y, a.z, b.z, c.z];
      const det = (
        M[0]*(M[4]*M[8]-M[5]*M[7]) -
        M[1]*(M[3]*M[8]-M[5]*M[6]) +
        M[2]*(M[3]*M[7]-M[4]*M[6])
      );
      if (Math.abs(det) < 1e-12) return;
      const invDet = 1/det;
      const inv = [
        (M[4]*M[8]-M[5]*M[7])*invDet,
        (M[2]*M[7]-M[1]*M[8])*invDet,
        (M[1]*M[5]-M[2]*M[4])*invDet,
        (M[5]*M[6]-M[3]*M[8])*invDet,
        (M[0]*M[8]-M[2]*M[6])*invDet,
        (M[2]*M[3]-M[0]*M[5])*invDet,
        (M[3]*M[7]-M[4]*M[6])*invDet,
        (M[1]*M[6]-M[0]*M[7])*invDet,
        (M[0]*M[4]-M[1]*M[3])*invDet
      ];
      function toFrac(pos) {
        const x=pos.x,y=pos.y,z=pos.z;
        return {
          u: inv[0]*x + inv[1]*y + inv[2]*z,
          v: inv[3]*x + inv[4]*y + inv[5]*z,
            w: inv[6]*x + inv[7]*y + inv[8]*z
        };
      }
      const fracs = molecule.atoms.map(a0 => ({ f: toFrac(a0.pos), atom: a0 }));
      // Apply new cell vectors
      if (typeof molecule.setCellVectors === 'function') {
        molecule.setCellVectors(newA, newB, newC);
      } else {
        cs.vectors.a = newA; cs.vectors.b = newB; cs.vectors.c = newC;
      }
      // Rebuild quick matrix for fractional to cartesian
      const A = newA, B = newB, C = newC;
      function fracToCart(f) { return {
        x: A.x*f.u + B.x*f.v + C.x*f.w,
        y: A.y*f.u + B.y*f.v + C.y*f.w,
        z: A.z*f.u + B.z*f.v + C.z*f.w
      }; }
      for (const entry of fracs) {
        const cart = fracToCart(entry.f);
        entry.atom.pos.x = cart.x;
        entry.atom.pos.y = cart.y;
        entry.atom.pos.z = cart.z;
      }
      if (molecule.changeCounter !== undefined) molecule.changeCounter++;
      updateCallback();
      const newVol = newA.x*(newB.y*newC.z - newB.z*newC.y) - newA.y*(newB.x*newC.z - newB.z*newC.x) + newA.z*(newB.x*newC.y - newB.y*newC.x);
      const dVpct = ((newVol - oldVol)/oldVol)*100;
      console.log('[Relaxation][Stress] applied cell update', {
        stressVoigt: stress.slice(0,6).map(v=> Number.isFinite(v)? Number(v.toFixed(5)) : v),
        hydrostatic_p: Number(p.toFixed(6)),
        scaleStep: Number(scaleStep.toFixed(6)),
        oldVolume: Number(oldVol.toFixed(6)),
        newVolume: Number(newVol.toFixed(6)),
        dV_percent: Number(dVpct.toFixed(4)),
        newLengths: {
          a: Number(newA.length().toFixed(4)),
          b: Number(newB.length().toFixed(4)),
          c: Number(newC.length().toFixed(4))
        }
      });
    } catch(e) {
      console.warn('[Relaxation] cell relax failed', e.message);
    }
  }

  async function performOptimizationStep() {
    const { energy, forces, stress } = await computeEnergyAndForces();
    
    // Check convergence
    const maxForce = getMaxForce(forces);
    const energyChange = lastEnergy !== null ? Math.abs(energy - lastEnergy) : Infinity;
    
    state.step = currentStep;
    state.energy = energy;
    state.maxForce = maxForce;
    state.energyChange = energyChange;

    console.log(`[Relaxation] Step ${currentStep}: E=${energy.toFixed(6)}, ΔE=${energyChange.toFixed(6)}, MaxF=${maxForce.toFixed(6)}`);

    // Check convergence criteria
    if (maxForce < forceTolerance && energyChange < energyTolerance) {
      state.converged = true;
      state.status = 'converged';
      isRunning = false;
      console.log(`[Relaxation] ✅ CONVERGED at step ${currentStep}! Energy change (${energyChange.toFixed(6)}) < tolerance (${energyTolerance})`);
      return true; // converged
    }

    // Compute search direction based on optimizer
    let direction;
    
    if (optimizer === 'steepest_descent') {
      // Simple steepest descent: direction = forces (negative gradient)
      direction = forces.map(f => vectorScale(f, 1.0));
    } else if (optimizer === 'conjugate_gradient') {
      if (searchDirection === null || lastForces === null) {
        // First step: use steepest descent
        direction = forces.map(f => vectorScale(f, 1.0));
      } else {
        // Compute Fletcher-Reeves beta
        const forcesDotForces = forces.reduce((sum, f) => sum + vectorDot(f, f), 0);
        const lastForcesDotLastForces = lastForces.reduce((sum, f) => sum + vectorDot(f, f), 0);
        beta = forcesDotForces / Math.max(lastForcesDotLastForces, 1e-12);
        
        // Compute new direction: forces + beta * old_direction
        direction = forces.map((f, i) => 
          vectorAdd(f, vectorScale(searchDirection[i], beta))
        );
      }
      searchDirection = direction;
    }

    // Normalize direction
    const directionMagnitude = Math.sqrt(direction.reduce((sum, d) => sum + vectorDot(d, d), 0));
    if (directionMagnitude < 1e-12) {
      state.status = 'converged';
      state.converged = true;
      isRunning = false;
      return true;
    }
    
    direction = direction.map(d => vectorScale(d, 1.0 / directionMagnitude));

    // Update positions
    const currentPositions = getAtomPositions();
    const newPositions = currentPositions.map((pos, i) => 
      vectorAdd(pos, vectorScale(direction[i], stepSize))
    );

    setAtomPositions(newPositions);
    // Post-move wrap to keep atoms inside cell
    try {
      if (molecule && molecule.wrapIntoCell) {
        const wrapped = molecule.wrapIntoCell({ refresh: false, bondRefresh: false, log: false });
        if (wrapped) updateCallback();
      }
    } catch {}

    // Record energy step for plotting (same interface as user interactions)
    try {
      energyCallback(energy);
    } catch (e) {
      // Ignore plot errors to not break relaxation
      console.warn('[Relaxation] Energy callback error:', e);
    }

    // Store for next iteration
    lastEnergy = energy;
    lastForces = forces;
    // Attempt simple cell relaxation using hydrostatic stress component
    if (stress && Array.isArray(stress) && stress.length >= 3) {
      tryCellRelax(stress);
    }
    currentStep++;

    return false; // not converged
  }

  async function start() {
    if (isRunning) return false;
    
    console.log(`[Relaxation] Starting relaxation with ${optimizer} optimizer`);
    reset();
    isRunning = true;
    stopRequested = false;
    currentStep = 0;
    state.status = 'running';
    
    try {
      while (isRunning && !stopRequested && currentStep < maxSteps) {
        const converged = await performOptimizationStep();
        if (converged) break;
        
        currentStep++;
        
        // Add minimal delay to allow UI updates and prevent blocking
        await new Promise(resolve => setTimeout(resolve, Math.max(1, stepDelay)));
      }
      
      if (stopRequested) {
        state.status = 'stopped';
      } else if (currentStep >= maxSteps && !state.converged) {
        state.status = 'max_steps';
        console.log(`[Relaxation] Reached maximum steps (${maxSteps})`);
      }
    } catch (error) {
      console.error('[Relaxation] Error during optimization:', error);
      state.status = 'error';
      state.error = error.message;
    } finally {
      isRunning = false;
    }
    
    return state.converged;
  }

  function stop() {
    if (isRunning) {
      console.log(`[Relaxation] Stop requested at step ${currentStep}`);
      stopRequested = true;
      state.status = 'stopped';
    }
    isRunning = false;
  }

  async function singleStep() {
    if (isRunning) return false;
    
    try {
      state.status = 'running';
      const converged = await performOptimizationStep();
      if (!converged && currentStep < maxSteps) {
        state.status = 'idle';
      } else if (currentStep >= maxSteps) {
        state.status = 'max_steps';
      }
      return converged;
    } catch (error) {
      console.error('[Relaxation] Error during single step:', error);
      state.status = 'error';
      state.error = error.message;
      return false;
    }
  }

  function getCurrentState() {
    return { ...state };
  }

  // Configuration setters
  function setStepSize(newStepSize) {
    stepSize = Math.max(0.001, Math.min(1.0, newStepSize));
  }

  function setMaxSteps(newMaxSteps) {
    maxSteps = Math.max(1, newMaxSteps);
  }

  function setOptimizer(newOptimizer) {
    if (['steepest_descent', 'conjugate_gradient'].includes(newOptimizer)) {
      optimizer = newOptimizer;
      // Reset CG state when switching optimizers
      searchDirection = null;
      beta = 0;
    }
  }

  function setTolerances(energy, force) {
    if (typeof energy === 'number' && energy > 0) {
      energyTolerance = energy;
    }
    if (typeof force === 'number' && force > 0) {
      forceTolerance = force;
    }
  }

  function setStepDelay(delay) {
    stepDelay = Math.max(0, delay);
  }

  return {
    start,
    stop,
    singleStep,
    reset,
    getCurrentState,
    setStepSize,
    setMaxSteps,
    setOptimizer,
    setTolerances,
    setStepDelay,
    get isRunning() { return isRunning; },
    get stepSize() { return stepSize; },
    get maxSteps() { return maxSteps; },
    get optimizer() { return optimizer; },
    get energyTolerance() { return energyTolerance; },
    get forceTolerance() { return forceTolerance; },
    get stepDelay() { return stepDelay; }
  };
}
