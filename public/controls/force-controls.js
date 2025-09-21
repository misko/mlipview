// controls/force-controls.js - Force visualization controls
export function createForceControls(forceVis) {
  let showForces = true; // Default to ON
  let lengthMode = "normalized";

  // Initialize forces as enabled
  forceVis.setEnabled(showForces);

  function toggleForces() {
    showForces = !showForces;
    forceVis.setEnabled(showForces);
    return showForces;
  }

  function toggleLengthMode() {
    lengthMode = (lengthMode === "normalized") ? "true" : "normalized";
    forceVis.setMode(lengthMode);
    return lengthMode;
  }

  function updateForces(forces) {
    if (showForces) {
      forceVis.setForces(forces);
    }
  }

  // Expose setForceScale globally for console access
  window.setForceScale = (s) => forceVis.setScaleTrue(s);

  return {
    isEnabled: () => showForces,
    getLengthMode: () => lengthMode,
    toggleForces,
    toggleLengthMode,
    updateForces
  };
}
