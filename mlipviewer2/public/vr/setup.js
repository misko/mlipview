// VR scaffolding placeholder: wraps Babylon WebXR experience when available.
// Provides a minimal interface so other modules can request entering/exiting VR.

export function createVRSupport(scene, { picking } = {}) {
  let xrHelper = null;
  let supported = false;
  function wireControllerEvents() {
    if (!xrHelper || !xrHelper.input || !picking) return;
    const input = xrHelper.input;
    // Listen for new controllers; map primary button/trigger to a pick in view direction.
    input.onControllerAddedObservable.add(controller => {
      const gp = controller.motionController;
      if (!gp) return;
      // primary component (trigger) might be named differently across profiles; scan for first with 'trigger' or 'select'
      const entries = Object.values(gp.components || {});
      const triggerComp = entries.find(c => /trigger|select/i.test(c.type || c.id || '')) || entries[0];
      if (!triggerComp) return;
      triggerComp.onButtonStateChangedObservable.add(comp => {
        if (comp.changes.pressed && comp.pressed) {
          // Perform a forward ray pick from controller
          const ray = controller.getForwardRay?.(100.0);
          let pickRes = null;
          if (ray && scene.pickWithRay) {
            pickRes = scene.pickWithRay(ray);
          }
          if (pickRes && pickRes.hit) {
            // Emulate pointer selection using existing picking logic: temporarily set pointer coords if needed
            const resolved = (function(resolvePick){
              const atom = resolvePick.resolveAtomPick(pickRes); if (atom) return atom;
              const bond = resolvePick.resolveBondPick(pickRes); if (bond) return bond;
              return null;
            })(picking.view || picking);
            if (resolved) {
              const selSvc = picking.selectionService || picking.selection || picking;
              if (resolved.kind === 'atom' && selSvc.clickAtom) selSvc.clickAtom(resolved.index);
              else if (resolved.kind === 'bond' && selSvc.clickBond) selSvc.clickBond(resolved);
              else if (selSvc.clear) selSvc.clear();
            } else if (picking.selectionService?.clear) {
              picking.selectionService.clear();
            }
          }
        }
      });
    });
  }
  async function init() {
    if (!scene.createDefaultXRExperienceAsync) {
      return { supported:false };
    }
    try {
      xrHelper = await scene.createDefaultXRExperienceAsync({});
      supported = true;
      wireControllerEvents();
    } catch (e) {
      console.warn('[vr] XR init failed', e);
      supported = false;
    }
    return { supported };
  }
  function isSupported() { return supported; }
  function enterVR() { if (xrHelper?.baseExperience) xrHelper.baseExperience.enterXR(); }
  function exitVR() { if (xrHelper?.baseExperience) xrHelper.baseExperience.exitXR(); }
  return { init, isSupported, enterVR, exitVR };
}
