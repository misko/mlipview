// mlipviewer2/public/interaction/adapters/vrInputAdapter.js
// Stub VR adapter: Provides interface parity; actual XR controller integration to be implemented.
// For now, we simulate periodic ray picks when enabled.

export function createVRInputAdapter({ scene, dragCore, pickingFacade, selectionService }) {
  let enabled = false;
  let lastPickTime = 0;
  const PICK_INTERVAL_MS = 150; // throttle simulated hover/pick

  function enable() { enabled = true; }
  function disable() { enabled = false; }

  // Simulated update loop using onBeforeRenderObservable
  if (!scene.onBeforeRenderObservable) {
    scene.onBeforeRenderObservable = { _c:[], add(fn){ this._c.push(fn); }, run(){ for(const f of this._c) try{f();}catch(_){} } };
  }
  scene.onBeforeRenderObservable.add(() => {
    if (!enabled) return;
    const now = performance.now ? performance.now() : Date.now();
    if (now - lastPickTime < PICK_INTERVAL_MS) return;
    lastPickTime = now;
    // Try bond first (mirroring desktop logic that prioritizes atom when clicked but here we just simulate):
    const bond = pickingFacade.pickBondAtPointer();
    if (bond) {
      selectionService.selectBond({ i: bond.i, j: bond.j, orientation: 0 });
      return;
    }
    const atom = pickingFacade.pickAtomAtPointer();
    if (atom) {
      selectionService.selectAtom(atom.index);
    }
  });

  return { enable, disable, isEnabled: ()=>enabled };
}
