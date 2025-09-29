// mlipviewer2/public/ui/rotationControllerCore.js
// Core bond rotation controller independent of UI implementation.

export function createRotationControllerCore({ selectionService, manipulationService, defaultStepDeg = 5 }) {
  let stepDeg = defaultStepDeg;
  function rotate(sign) {
    const sel = selectionService.getSelection();
    if (!sel || sel.kind !== 'bond') return false;
    const { i, j, orientation } = sel.data;
    manipulationService.rotateBond({ i, j, orientation, angleDeg: sign * stepDeg });
    return true;
  }
  function setStep(v) { if (v > 0) stepDeg = v; }
  function getStep() { return stepDeg; }
  function cycleOrientation() {
    const sel = selectionService.getSelection();
    if (!sel || sel.kind !== 'bond') return;
    const { i, j, orientation } = sel.data;
    const next = (orientation === 0) ? 1 : 0;
    selectionService.selectBond({ i, j, orientation: next });
  }
  return { rotateMinus: ()=>rotate(-1), rotatePlus: ()=>rotate(1), setStep, getStep, cycleOrientation };
}
