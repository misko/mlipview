import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('x-vr trigger selection', () => {
  let currentTime;

  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
    currentTime = 0;
    global.performance = { now: () => currentTime };
  });

  afterEach(() => {
    delete global.window;
    delete global.document;
    delete global.navigator;
  });

  test('controller trigger picks atom via resolveAtomPick fallback', async () => {
    const controllerAdded = [];
    const controllers = [];

    const triggerComponent = { pressed: false };
    const pointerNode = {
      getAbsolutePosition: () => new BABYLON.Vector3(0, 0, 0),
      getDirection: () => new BABYLON.Vector3(0, 0, -1),
    };
    const controller = {
      uniqueId: 'ctrl-1',
      pointer: pointerNode,
      motionController: {
        getComponent(name) {
          return name === 'xr-standard-trigger' ? triggerComponent : null;
        },
      },
      inputSource: { gamepad: { buttons: [{ pressed: false }] } },
    };

    const xrHelper = {
      input: {
        controllers,
        onControllerAddedObservable: {
          add(fn) {
            controllerAdded.push(fn);
          },
        },
        onControllerRemovedObservable: { add() {} },
      },
      baseExperience: {
        sessionManager: {
          onXRSessionInit: { add() {} },
          onXRSessionEnded: { add() {} },
        },
        camera: {
          position: new BABYLON.Vector3(0, 1.6, 0),
          getDirection: () => new BABYLON.Vector3(0, 0, 1),
        },
      },
      pointerSelection: {},
      featuresManager: {
        enableFeature: jest.fn(),
      },
    };

    const engine = {
      setHardwareScalingLevel() {},
      getRenderingCanvas() {
        return { style: {} };
      },
    };

    const frameCallbacks = [];
    const scene = {
      uid: 7,
      meshes: [],
      onBeforeRenderObservable: {
        add(fn) {
          frameCallbacks.push(fn);
        },
      },
      getEngine() {
        return engine;
      },
      pickWithRay: jest.fn(() => ({
        hit: true,
        distance: 0.4,
        pickedMesh: { name: 'atom_5' },
      })),
    };
    scene.activeCamera = {
      position: new BABYLON.Vector3(0, 1.6, 0),
      getDirection: () => new BABYLON.Vector3(1, 0, 0),
    };

    global.window = { addEventListener() {} };
    global.document = {
      getElementById: () => ({ style: {} }),
      body: { style: {}, appendChild() {} },
    };
    global.navigator = { xr: { isSessionSupported: async () => false } };

    jest.unstable_mockModule('../public/vr/vr-utils.js', () => ({
      getMoleculeMasters: jest.fn(() => []),
      refreshMoleculeMasters: jest.fn(() => []),
      getMoleculeDiag: jest.fn(() => 2),
      getAnyMaster: jest.fn(() => null),
      transformLocalToWorld: jest.fn((_scene, vec) => ({ ...vec })),
    }));
    jest.unstable_mockModule('../public/vr/vr-laser.js', () => ({
      createVRLaserManager: jest.fn(() => ({
        ensureForController: jest.fn(),
        updateFrame: jest.fn(),
        getInfo: jest.fn(() => []),
      })),
    }));
    jest.unstable_mockModule('../public/vr/vr-picker.js', () => ({
      createVRPicker: jest.fn(() => ({
        pickAtomWithRay: () => ({ idx: 5 }),
        pickBondWithRay: () => null,
      })),
    }));
    jest.unstable_mockModule('../public/vr/xr-hud-bars.js', () => ({
      ensureWorldHUD: jest.fn(() => false),
      forceWorldHUD: jest.fn(() => false),
    }));

    scene.createDefaultXRExperienceAsync = jest.fn(async () => xrHelper);

    const selectionService = { clickAtom: jest.fn(), clickBond: jest.fn() };
    const view = {
      resolveAtomPick: jest.fn(() => ({ idx: 5 })),
      resolveBondPick: jest.fn(() => null),
    };
    const picking = { view, selectionService };

    const { createVRSupport } = await import('../public/vr/setup.js');
    const vr = createVRSupport(scene, { picking });
    await vr.init();

    // Simulate controller connection.
    controllers.push(controller);
    controllerAdded.forEach((fn) => fn(controller));

    expect(frameCallbacks.length).toBeGreaterThan(0);
    const tick = frameCallbacks[0];

    // Frame with trigger pressed to drive selection.
    currentTime = 20;
    tick();
    const beforeCalls = scene.pickWithRay.mock.calls.length;

    triggerComponent.pressed = true;
    controller.inputSource.gamepad.buttons[0].pressed = true;
    currentTime = 50;
    tick();

    const afterCalls = scene.pickWithRay.mock.calls.length;
    expect(afterCalls).toBeGreaterThan(beforeCalls);
  });
});
