// mlipviewer2/public/interaction/pickingFacade.js
// Wrap scene picking mechanics so adapters can uniformly request atom/bond picks.

export function createPickingFacade({ scene, moleculeView }) {
  // Expect moleculeView provides resolveAtomPick & resolveBondPick helpers.
  function pickAtomAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY);
    const res = moleculeView.resolveAtomPick(pick);
    return res ? { index: res.index, element: res.element } : null;
  }
  function pickBondAtPointer() {
    const pick = scene.pick(scene.pointerX, scene.pointerY);
    const res = moleculeView.resolveBondPick(pick);
    return res ? { i: res.i, j: res.j, key: res.key, index: res.index } : null;
  }
  return { pickAtomAtPointer, pickBondAtPointer };
}
