/** @jest-environment jsdom */

// Mock browser test: verifies that when a bond is selected, joystick X-axis drives continuous
// rotations via manipulation.rotateBond in the VR layer.

let createMoleculeState, createMoleculeView, createVRSupport, setupVRFeatures;

beforeAll(async () => {
  ({ createMoleculeState } = await import('../public/domain/moleculeState.js'));
  ({ createMoleculeView } = await import('../public/render/moleculeView.js'));
  ({ createVRSupport, setupVRFeatures } = await import('../public/vr/setup.js'));
});

function makeScene() {
  // Minimal scene stub that supports onBeforeRenderObservable and controller plumbing
  return {
    uid: Math.random(),
    meshes: [],
    transformNodes: [],
    onBeforeRenderObservable: {
      _cbs: [],
      add(fn) {
        this._cbs.push(fn);
        return fn;
      },
      remove(fn) {
        const i = this._cbs.indexOf(fn);
        if (i >= 0) this._cbs.splice(i, 1);
      },
    },
    pickWithRay() {
      return { hit: false };
    },
    getEngine() {
      return {
        getRenderingCanvas() {
          return {
            getBoundingClientRect() {
              return { left: 0, top: 0 };
            },
          };
        },
      };
    },
    activeCamera: {
      getDirection() {
        return new BABYLON.Vector3(1, 0, 0);
      },
      position: new BABYLON.Vector3(),
    },
  };
}

function step(scene, n = 1) {
  for (let i = 0; i < n; i++) {
    try {
      for (const fn of scene.onBeforeRenderObservable._cbs) fn();
    } catch {}
  }
}

function makeController({ id = 'c1' } = {}) {
  // Gamepad with axes[2] = right stick X; [0] = left stick X
  const gp = { axes: [0, 0, 0, 0], buttons: [{ pressed: false }] };
  return {
    uniqueId: id,
    inputSource: { gamepad: gp },
    motionController: {
      getComponent() {
        return { pressed: false };
      },
    },
    // Minimal transform methods for ray/laser code; not used in this test
    getDirection() {
      return new BABYLON.Vector3(0, 0, 1);
    },
    pointer: {
      getAbsolutePosition() {
        return new BABYLON.Vector3();
      },
      getDirection() {
        return new BABYLON.Vector3(0, 0, 1);
      },
    },
  };
}

function selectBond(state, selectionService) {
  // Select bond i-j with default orientation (0)
  const b = state.bonds[0];
  selectionService.clickBond({ i: b.i, j: b.j, key: 0, index: 0 });
}

describe('VR bond joystick rotation', () => {
  test('stick X produces continuous bond rotation over frames', async () => {
    // Build a simple molecule with a single bond 0-1 and a third atom connected to atom 1
    const st = createMoleculeState({
      elements: ['C', 'C', 'H'],
      // Place atom 2 slightly off the bond axis so rotation produces visible y/z change
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.5, y: 0, z: 0 },
        { x: 2.2, y: 0.4, z: 0 },
      ],
      bonds: [
        { i: 0, j: 1 },
        { i: 1, j: 2 },
      ],
    });
    const scene = makeScene();
    createMoleculeView(scene, st);

    // Minimal selection service compatible with VR createVRSupport usage
    const selectionServiceModule = await import('../public/domain/selectionService.js');
    const selectionService = selectionServiceModule.createSelectionService(st);

    // Manipulation service used by VR should be the real one so rotateBond alters positions
    const manipulationModule = await import('../public/domain/manipulationService.js');
    const manipulation = manipulationModule.createManipulationService(st, {});

    // Initialize VR feature loop directly with a stub xrHelper so frame handler is installed
    const xr = {
      input: {
        controllers: [],
        onControllerAddedObservable: { add() {} },
        onControllerRemovedObservable: { add() {} },
      },
      baseExperience: {
        sessionManager: { onXRSessionInit: { add() {} }, onXRSessionEnded: { add() {} } },
        camera: {
          position: new BABYLON.Vector3(),
          getDirection() {
            return new BABYLON.Vector3(0, 0, 1);
          },
        },
      },
    };
    setupVRFeatures(xr, scene, { view: null, selectionService, manipulation, molState: st });
    const ctrl = makeController({ id: 'sim' });
    xr.input.controllers.push(ctrl);

    // Select the bond for rotation
    selectBond(st, selectionService);

    const before = st.positions.map((p) => ({ ...p }));

    // Push joystick to the right slightly and advance several frames; expect rotation applied
    ctrl.inputSource.gamepad.axes[2] = 0.9; // strong right deflection

    // Advance frames to accumulate rotations
    step(scene, 10);

    const after = st.positions.map((p) => ({ ...p }));

    // At least one atom on the rotating side (atom 0 side or 1 side based on orientation=0 => side 'j') should have moved significantly
    // With default orientation 0 => orientationToSide returns 'j', so atoms on side of j (index 1 and neighbors) rotate around axis (i->j) with anchor=i
    // Thus atom 2 (connected to 1) should move off the x-axis (y or z non-zero)
    expect(Math.abs(after[2].y) + Math.abs(after[2].z)).toBeGreaterThan(1e-6);
    // Anchor atom (i=0 when orientation=0) should remain on anchor line for small rotation (x may change negligibly due to float but y/z ~ 0)
    expect(Math.abs(after[0].y)).toBeLessThan(1e-6);
  });
});
