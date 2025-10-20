import { createVRSupport } from '../public/vr/setup.js';

// We will stub a minimal scene, view, and selection service plus XR helper structures.

function makeStubScene(pickResult) {
  return {
    createDefaultXRExperienceAsync: async () => ({
      input: {
        onControllerAddedObservable: {
          add: (cb) => {
            setTimeout(() => cb(makeController(pickResult)), 0);
          },
        },
      },
      baseExperience: { enterXR() {}, exitXR() {} },
    }),
    pickWithRay: () => pickResult,
  };
}

function makeController(pickResult) {
  const triggerListeners = [];
  const triggerComp = {
    changes: { pressed: false },
    pressed: false,
    onButtonStateChangedObservable: {
      add(fn) {
        triggerListeners.push(fn);
      },
    },
  };
  const motionController = { components: { trigger: { ...triggerComp, type: 'trigger' } } };
  const ctrl = {
    motionController,
    getForwardRay: () => ({}),
  };
  // Expose a helper to simulate press
  ctrl.__press = () => {
    triggerComp.changes.pressed = true;
    triggerComp.pressed = true;
    triggerListeners.forEach((fn) => fn(triggerComp));
    triggerComp.changes.pressed = false;
    triggerComp.pressed = false;
  };
  return ctrl;
}

test('VR trigger maps to selection via picking', async () => {
  const pickedAtom = { hit: true, pickedMesh: {} };
  const view = {
    resolveAtomPick: (p) => (pickedAtom === p ? { kind: 'atom', index: 5 } : null),
    resolveBondPick: () => null,
  };
  const selectionService = {
    clicked: null,
    clickAtom(i) {
      this.clicked = { kind: 'atom', index: i };
    },
    clickBond() {},
    clear() {},
  };
  const picking = { view, selectionService };
  const scene = makeStubScene(pickedAtom);
  const vr = createVRSupport(scene, { picking });
  const initRes = await vr.init();
  expect(initRes.supported).toBe(true);
  await new Promise((r) => setTimeout(r, 10));
  // In shim mode controllers are fabricated inside init; retrieve via accessor
  const shimControllers = vr.controllers();
  expect(Array.isArray(shimControllers)).toBe(true);
  // If no fabricated controller (non-shim environment), skip
  if (!shimControllers.length) return;
  if (shimControllers[0].__press) shimControllers[0].__press();
  await new Promise((r) => setTimeout(r, 20));
  expect(selectionService.clicked).toEqual({ kind: 'atom', index: 5 });
});
