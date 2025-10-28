import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('x-vr controller handshake', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
  });

  afterEach(() => {
    delete global.window;
    delete global.document;
    delete global.navigator;
    delete global.performance;
  });

  test('controllers are surfaced and fed into laser manager each frame', async () => {
    const controllers = [];
    const controllerAdded = [];

    const pointer = {
      getAbsolutePosition: () => new BABYLON.Vector3(0, 0, 0),
      getDirection: () => new BABYLON.Vector3(0, 0, -1),
    };
    const controller = {
      uniqueId: 'ctrl-handshake',
      pointer,
      motionController: { getComponent: () => ({ pressed: false }) },
      inputSource: { gamepad: { buttons: [{ pressed: false }], axes: [] } },
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
      },
      pointerSelection: {},
      featuresManager: { enableFeature: jest.fn() },
    };

    const engine = {
      setHardwareScalingLevel() {},
      getRenderingCanvas() {
        return { style: {} };
      },
    };
    const frameCallbacks = [];
    const scene = {
      uid: 303,
      meshes: [],
      onBeforeRenderObservable: {
        add(fn) {
          frameCallbacks.push(fn);
        },
      },
      getEngine() {
        return engine;
      },
      pickWithRay: () => ({ hit: false }),
    };

    global.performance = { now: () => 0 };
    global.window = { addEventListener() {} };
    global.document = {
      getElementById: () => ({ style: {} }),
      body: { style: {} },
    };
    global.navigator = { xr: { isSessionSupported: async () => false } };

    await jest.unstable_mockModule('../public/vr/vr-utils.js', () => ({
      getMoleculeMasters: jest.fn(() => []),
      refreshMoleculeMasters: jest.fn(() => []),
      getMoleculeDiag: jest.fn(() => 1),
      getAnyMaster: jest.fn(() => null),
      transformLocalToWorld: jest.fn((_scene, v) => ({ ...v })),
    }));
    const createVRLaserManagerMock = jest.fn(() => ({
      ensureForController: jest.fn(),
      updateFrame: jest.fn(),
      getInfo: jest.fn(() => []),
    }));
    await jest.unstable_mockModule('../public/vr/vr-laser.js', () => ({
      createVRLaserManager: createVRLaserManagerMock,
    }));
    await jest.unstable_mockModule('../public/vr/vr-picker.js', () => ({
      createVRPicker: jest.fn(() => null),
    }));
    await jest.unstable_mockModule('../public/vr/xr-hud-bars.js', () => ({
      ensureWorldHUD: jest.fn(() => false),
      forceWorldHUD: jest.fn(() => false),
    }));

    scene.createDefaultXRExperienceAsync = jest.fn(async () => xrHelper);

    let createVRSupport;
    await jest.isolateModulesAsync(async () => {
      ({ createVRSupport } = await import('../public/vr/setup.js'));
    });
    const vr = createVRSupport(scene, { picking: {} });
    await vr.init();

    controllers.push(controller);
    controllerAdded.forEach((fn) => fn(controller));

    expect(frameCallbacks.length).toBeGreaterThan(0);
    frameCallbacks[0]();

    expect(vr.controllers()).toEqual([controller]);
  });
});
