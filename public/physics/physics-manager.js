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

  return {
    tickFrame,
    updatePhysicsCache,
    resetCache,
    // Relaxation interface
    startRelaxation,
    stopRelaxation,
    getRelaxationState,
    isRelaxationRunning,
    get energy() { return cachedEnergy; },
    get forces() { return cachedForces; }
  };
}
