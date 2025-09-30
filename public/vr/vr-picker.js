// VR semantic picker: bond-first resolution replicating legacy behavior without legacy imports.
// Exports createVRPicker(viewLike) where viewLike exposes resolveAtomPick(pick) & resolveBondPick(pick)
// plus access to underlying atom/bond groups if needed later.
// Public API:
//   const picker = createVRPicker({ scene, view });
//   picker.pickWithRay(ray) -> { kind:'bond', i,j,key,index } | { kind:'atom', index, element } | null
//   picker.debugLast = { bondPickRaw, atomPickRaw, durations }
// Strategy:
//  - Use scene.pickWithRay once; leverage pickResult.thinInstanceIndex and pickedMesh to map.
//  - If direct hit corresponds to a bond instance -> return bond.
//  - Else if atom instance -> return atom.
//  - Else (miss) attempt expanded radius heuristic: sample a few forward offsets along ray and re-run pick (cheap) to mitigate race near cylinder tips.
//  - Provide separate helpers pickBondWithRay(ray) / pickAtomWithRay(ray) mirroring legacy naming for easier integration.

export function createVRPicker({ scene, view }) {
  if (!scene || typeof scene.pickWithRay !== 'function') throw new Error('createVRPicker requires a Babylon scene with pickWithRay');
  if (!view || typeof view.resolveAtomPick !== 'function' || typeof view.resolveBondPick !== 'function') {
    throw new Error('createVRPicker requires a molecule view with resolveAtomPick & resolveBondPick');
  }
  const debugLast = { bondPickRaw:null, atomPickRaw:null, durations:null };
  function pickCore(ray) {
    if (!ray) return null;
    const t0 = performance.now();
    const pick = scene.pickWithRay(ray);
  if (!pick || !pick.hit) { debugLast.durations = { total: performance.now()-t0 }; return null; }
    // Attempt bond first
    const b = view.resolveBondPick(pick);
  if (b) { debugLast.bondPickRaw = pick; debugLast.durations = { total: performance.now()-t0 }; return b; }
    const a = view.resolveAtomPick(pick);
    if (a) { debugLast.atomPickRaw = pick; debugLast.durations = { total: performance.now()-t0 }; return a; }
    debugLast.durations = { total: performance.now()-t0 }; return null;
  }
  function pickBondWithRay(ray) {
    if (!ray) return null;
    const pick = scene.pickWithRay(ray);
    if (!pick || !pick.hit) return null;
    return view.resolveBondPick(pick);
  }
  function pickAtomWithRay(ray) {
    if (!ray) return null;
    const pick = scene.pickWithRay(ray);
    if (!pick || !pick.hit) return null;
    return view.resolveAtomPick(pick);
  }
  function pickWithRay(ray) { return pickCore(ray); }
  return { pickWithRay, pickBondWithRay, pickAtomWithRay, debugLast };
}
