// Unified physics manager for desktop & VR to keep logic in sync.
// Provides cached compute, throttled retries, and backend switching reset.

import { createRelaxationEngine } from './relaxation.js';

export function createPhysicsManager({ state, getMlip, energyPlot }) {
  let cachedEnergy = null;
  let cachedForces = null;
  let lastChangeCounter = -1;
  let inFlight = false;
  let frameCount = 0;
  let lastComputeAttemptFrame = -Infinity;
  let computeSeq = 0;

  const MIN_RETRY_INTERVAL = 30; // frames before retry when energy null

  // Relaxation engine
  let relaxationEngine = null;
  // MD engine state
  let mdRunning = false;
  let mdStep = 0;
  let mdVelocities = null; // Array of BABYLON.Vector3
  let mdHandle = null; // requestAnimationFrame id
  let mdOpts = null;

  const MASS_TABLE = { H: 1.0079, C: 12.011, N: 14.007, O: 15.999, S: 32.06 };
  function atomMass(el) { return MASS_TABLE[el] || 12.0; }
  
  function ensureRelaxationEngine() {
    if (!relaxationEngine) {
      relaxationEngine = createRelaxationEngine({
        getMlip,
        molecule: state.molecule,
        updateCallback: () => {
          // Update visual representation - atoms and chemistry-aware bonds
          if (state.molecule && typeof state.molecule.refreshAtoms === 'function') {
            state.molecule.refreshAtoms();
          }
          if (state.molecule && typeof state.molecule.refreshBonds === 'function') {
            state.molecule.refreshBonds(); // Now includes chemistry-aware updating
          }
          
          // Trigger physics cache update when coordinates change during relaxation
          resetCache();
          updatePhysicsCache();
        },
        energyCallback: (energy) => {
          // Record energy step for plotting (same interface as user interactions)
          // Only available in desktop mode, VR doesn't have energy plotting
          if (energyPlot && 
              typeof energyPlot.isEnabled === 'function' && energyPlot.isEnabled() && 
              typeof energyPlot.recordStep === 'function') {
            energyPlot.recordStep(energy);
          }
        }
      });
    }
    return relaxationEngine;
  }

  function hasStructureChanged() {
    const molecule = state.molecule;
    if (molecule && molecule.changeCounter !== lastChangeCounter) {
      lastChangeCounter = molecule.changeCounter;
      return true;
    }
    return false;
  }

  function updatePhysicsCache() {
    const structureChanged = hasStructureChanged();
    if (structureChanged) {
      cachedEnergy = null;
      cachedForces = null;
    }
    if (inFlight) return false;
    const shouldAttempt = structureChanged || (cachedEnergy === null && (frameCount - lastComputeAttemptFrame) >= MIN_RETRY_INTERVAL);
    if (!shouldAttempt) return false;
    inFlight = true;
    lastComputeAttemptFrame = frameCount;
    const seq = ++computeSeq;
    let result;
    let mlip;
    try { mlip = getMlip(); result = mlip.compute(); } catch (e) { inFlight = false; return false; }
    if (result && typeof result.then === 'function') {
      result.then(res => {
        inFlight = false;
        if (seq !== computeSeq) return;
        if (res && Number.isFinite(res.energy)) cachedEnergy = res.energy; else cachedEnergy = null;
        cachedForces = Array.isArray(res?.forces) ? res.forces : null;
      }).catch(() => { inFlight = false; cachedEnergy = null; });
      return false;
    }
    inFlight = false;
    if (result) {
      if (Number.isFinite(result.energy)) cachedEnergy = result.energy; else cachedEnergy = null;
      cachedForces = Array.isArray(result.forces) ? result.forces : null;
      return true;
    }
    return false;
  }

  function tickFrame() { frameCount++; }

  function resetCache() {
    cachedEnergy = null;
    cachedForces = null;
    lastChangeCounter = -1;
    computeSeq++; // invalidate async
  }

  // Relaxation methods
  function startRelaxation(options = {}) {
    // Mutual exclusion: stop MD if running
    if (mdRunning) stopMD();
    const engine = ensureRelaxationEngine();
    // Update engine parameters if provided
    if (options.optimizer) engine.setOptimizer(options.optimizer);
    if (options.stepSize) engine.setStepSize(options.stepSize);
    if (options.maxSteps) engine.setMaxSteps(options.maxSteps);
    if (options.stepDelay !== undefined) engine.setStepDelay(options.stepDelay);
    if (options.energyTolerance || options.forceTolerance) {
      engine.setTolerances(options.energyTolerance, options.forceTolerance);
    }
    return engine.start();
  }

  function stopRelaxation() {
    if (relaxationEngine) {
      relaxationEngine.stop();
    }
  }



  function getRelaxationState() {
    if (!relaxationEngine) {
      return { status: 'idle', step: 0, energy: null, maxForce: null, converged: false };
    }
    return relaxationEngine.getCurrentState();
  }

  function isRelaxationRunning() {
    return relaxationEngine ? relaxationEngine.isRunning : false;
  }

  // --------------------------------------------------------------
  // MD (Molecular Dynamics) basic velocity Verlet integrator
  // --------------------------------------------------------------
  function initVelocities() {
    const atoms = state.molecule?.atoms || [];
    mdVelocities = atoms.map(() => new BABYLON.Vector3(0,0,0));
  }
  function isMDRunning() { return mdRunning; }
  function getMDState() { return { running: mdRunning, step: mdStep }; }

  function currentKineticTemperature() {
    const atoms = state.molecule?.atoms || [];
    if (!atoms.length || !mdVelocities) return 0;
    // Simplified k_B = 1 units for visualization: T = (2 KE)/(3N)
    let ke = 0;
    for (let i=0;i<atoms.length;i++) {
      const m = atomMass(atoms[i].element);
      const v = mdVelocities[i];
      ke += 0.5 * m * (v.x*v.x + v.y*v.y + v.z*v.z);
    }
    const dof = Math.max(1, 3*atoms.length - 3); // subtract rigid translation
    return (2*ke) / dof;
  }

  async function mdSingleStep() {
    if (!mdRunning) return;
    const mol = state.molecule;
    if (!mol) { stopMD(); return; }
    const atoms = mol.atoms;
    if (!atoms || !atoms.length) { stopMD(); return; }
    // Ensure current coordinates are wrapped before computing forces
    if (mol.wrapIntoCell) mol.wrapIntoCell({ refresh: true, bondRefresh: false, log: false });
    const mlip = getMlip();
    try {
      // 1. Compute forces at current positions
      const res = await mlip.compute();
      if (!res || !Array.isArray(res.forces)) throw new Error('MD: invalid forces');
      const forces = res.forces; // array of {x,y,z}
      const dt = mdOpts.dt;
      const half = dt * 0.5;
      // Velocity Verlet: v(t+1/2)
      for (let i=0;i<atoms.length;i++) {
        const m = atomMass(atoms[i].element);
        const f = forces[i];
        const v = mdVelocities[i];
        v.x += (f.x / m) * half;
        v.y += (f.y / m) * half;
        v.z += (f.z / m) * half;
      }
      // Update positions x(t+1)
      for (let i=0;i<atoms.length;i++) {
        const a = atoms[i];
        const v = mdVelocities[i];
        a.pos.x += v.x * dt;
        a.pos.y += v.y * dt;
        a.pos.z += v.z * dt;
      }
      // Wrap after drift
      if (mol.wrapIntoCell) mol.wrapIntoCell({ refresh: false, bondRefresh: false, log: false });
      if (mol.changeCounter !== undefined) mol.changeCounter++;
      if (typeof mol.refreshAtoms === 'function') mol.refreshAtoms();
      // Recompute bonds each step (connectivity may change if system distorts/explodes)
      try {
        if (typeof mol.recomputeBonds === 'function') {
          mol.recomputeBonds({ hysteresis: 1.01 }); // slightly tighter to avoid ghosty long bonds
        } else if (typeof mol.refreshBonds === 'function') {
          mol.refreshBonds();
        }
      } catch(e) { console.warn('[MD] bond update failed', e.message); }
      resetCache(); // positions changed
      // 2. Compute forces at new positions for full step velocity update
      if (mol.wrapIntoCell) mol.wrapIntoCell({ refresh: false, bondRefresh: false, log: false });
      const res2 = await mlip.compute();
      if (!res2 || !Array.isArray(res2.forces)) throw new Error('MD: invalid forces 2');
      const forces2 = res2.forces;
      for (let i=0;i<atoms.length;i++) {
        const m = atomMass(atoms[i].element);
        const f = forces2[i];
        const v = mdVelocities[i];
        v.x += (f.x / m) * half;
        v.y += (f.y / m) * half;
        v.z += (f.z / m) * half;
      }
      // --- Thermostat (Berendsen style) ---
      if (mdOpts.thermostat) {
        const Tcur = currentKineticTemperature();
        const Ttarget = mdOpts.thermostat.target;
        const tau = mdOpts.thermostat.tau; // relaxation time in MD time units
        if (Tcur > 0 && Ttarget > 0 && tau > 0) {
          const dt = mdOpts.dt;
          // Berendsen scaling factor: lambda = sqrt(1 + (dt/tau)*(Ttarget/Tcur - 1))
          const lambda2 = 1 + (dt / tau) * (Ttarget / Tcur - 1);
          if (lambda2 > 0) {
            const lambda = Math.sqrt(lambda2);
            for (const v of mdVelocities) {
              v.x *= lambda; v.y *= lambda; v.z *= lambda;
            }
          }
        }
      }
      mdStep++;
      if (mdOpts.energyPlot && Number.isFinite(res2.energy)) {
        try { mdOpts.energyPlot.recordStep(res2.energy); } catch {}
      }
      if (mdRunning) mdHandle = requestAnimationFrame(mdSingleStep);
    } catch (e) {
      console.warn('[MD] step error, stopping:', e.message);
      stopMD();
    }
  }

  function startMD(options = {}) {
    if (isRelaxationRunning()) stopRelaxation();
    if (mdRunning) return;
    mdOpts = { 
      dt: options.dt || 0.5,
      energyPlot: options.energyPlot || energyPlot,
      // Default thermostat ON unless explicitly disabled
      thermostat: (options.thermostat === false) ? { target: 0, tau: 0 } : (options.thermostat || { target: 0.5, tau: 50 })
    };
    initVelocities();
    mdRunning = true;
    mdStep = 0;
    mdHandle = requestAnimationFrame(mdSingleStep);
    console.log('[MD] started with dt=', mdOpts.dt);
  }

  function setMDThermostat({ target, tau }) {
    if (!mdOpts) return;
    if (target <= 0 || tau <= 0) {
      mdOpts.thermostat = { target: 0, tau: 0 }; // disable
    } else {
      mdOpts.thermostat = { target, tau };
    }
    console.log('[MD] thermostat set to', mdOpts.thermostat);
  }
  function stopMD() {
    if (!mdRunning) return;
    mdRunning = false;
    if (mdHandle) cancelAnimationFrame(mdHandle);
    mdHandle = null;
    console.log('[MD] stopped at step', mdStep);
  }

  return {
    tickFrame,
    updatePhysicsCache,
    resetCache,
    // Relaxation interface
    startRelaxation,
    stopRelaxation,
    getRelaxationState,
    isRelaxationRunning,
    // MD interface
    startMD,
    stopMD,
    isMDRunning,
    getMDState,
    setMDThermostat,
    get energy() { return cachedEnergy; },
    get forces() { return cachedForces; }
  };
}
