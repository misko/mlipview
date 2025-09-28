import { createVRSupport } from '../public/vr/setup.js';

// We will stub a minimal scene, view, and selection service plus XR helper structures.

function makeStubScene(pickResult) {
  return {
    createDefaultXRExperienceAsync: async () => ({
      input: {
        onControllerAddedObservable: { add: (cb) => { setTimeout(()=>cb(makeController(pickResult)),0); } }
      },
      baseExperience: { enterXR() {}, exitXR() {} }
    }),
    pickWithRay: () => pickResult,
  };
}

function makeController(pickResult) {
  const triggerListeners = [];
  const triggerComp = {
    changes:{ pressed:false },
    pressed:false,
    onButtonStateChangedObservable:{ add(fn){ triggerListeners.push(fn);} }
  };
  const motionController = { components: { trigger: { ...triggerComp, type:'trigger' } } };
  const ctrl = {
    motionController,
    getForwardRay: () => ({})
  };
  // Expose a helper to simulate press
  ctrl.__press = () => {
    triggerComp.changes.pressed = true; triggerComp.pressed = true;
    triggerListeners.forEach(fn=>fn(triggerComp));
    triggerComp.changes.pressed = false; triggerComp.pressed = false;
  };
  return ctrl;
}

test('VR trigger maps to selection via picking', async () => {
  const pickedAtom = { hit:true, pickedMesh:{}, /* custom properties consumed by resolveAtomPick below */ };
  const view = {
    resolveAtomPick: (p)=> pickedAtom===p? { kind:'atom', index:5 } : null,
    resolveBondPick: ()=> null
  };
  const selectionService = { clicked:null, clickAtom(i){ this.clicked={kind:'atom', index:i}; }, clickBond(){}, clear(){} };
  const picking = { view, selectionService };
  const scene = makeStubScene(pickedAtom);
  const vr = createVRSupport(scene, { picking });
  const initRes = await vr.init();
  expect(initRes.supported).toBe(true);
  // Controller added async; wait a tick then simulate press
  await new Promise(r=>setTimeout(r,5));
  // Find controller reference and simulate press
  // Our stub stores controller only inside the added callback; expose through closure hack by reinitializing? Instead, rebuild to capture.
  // Simplify: rebuild a controller manually and invoke listeners again.
  // For test brevity, we just assert selection after a manual simulation path.
  // Re-add a controller to capture its instance for pressing.
  const controllers = [];
  scene.createDefaultXRExperienceAsync = async () => ({ input:{ onControllerAddedObservable:{ add:(cb)=>{ const c=makeController(pickedAtom); controllers.push(c); cb(c);} } }, baseExperience:{ enterXR(){}, exitXR(){} } });
  // Re-init to attach new handler capturing controller.
  await vr.init();
  await new Promise(r=>setTimeout(r,5));
  controllers[0].__press();
  expect(selectionService.clicked).toEqual({ kind:'atom', index:5 });
});
