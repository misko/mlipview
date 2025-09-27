// Unified physics manager for desktop & VR to keep logic in sync.
// Provides cached compute, throttled retries, and backend switching reset.

export function createPhysicsManager({ state, getMlip }) {
  let cachedEnergy = null;
  let cachedForces = null;
  let lastChangeCounter = -1;
  let inFlight = false;
  let frameCount = 0;
  let lastComputeAttemptFrame = -Infinity;
  let computeSeq = 0;

  const MIN_RETRY_INTERVAL = 30; // frames before retry when energy null

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

  return {
    tickFrame,
    updatePhysicsCache,
    resetCache,
    get energy() { return cachedEnergy; },
    get forces() { return cachedForces; }
  };
}
