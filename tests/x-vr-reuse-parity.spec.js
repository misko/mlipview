import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('x-vr reuse parity', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
  });

  afterEach(() => {
    delete global.window;
    delete global.document;
    delete global.navigator;
  });

  test('initVRApp reuses existing viewer scene and engine', async () => {
    const scene = {
      meshes: Array.from({ length: 6 }, (_, i) => ({ name: 'atom_' + i })),
      createDefaultXRExperienceAsync: jest.fn(async () => ({
        baseExperience: { sessionManager: {} },
        input: { controllers: [], onControllerAddedObservable: { add() {} } },
      })),
    };
    const engine = {
      runRenderLoop: jest.fn(),
      getRenderingCanvas: () => ({ style: {} }),
      resize: jest.fn(),
    };

    global.window = {
      _viewer: { engine, scene },
      addEventListener() {},
    };
    global.document = {
      getElementById: () => null,
      querySelector: () => null,
      createElement: () => ({ style: {} }),
      body: { appendChild() {}, style: {} },
    };
    global.navigator = { xr: { isSessionSupported: async () => false } };

    const { initVRApp } = await import('../public/vr/main-vr.js');
    const res = await initVRApp();
    expect(res.scene).toBe(scene);
    expect(res.engine).toBe(engine);
    expect(scene.createDefaultXRExperienceAsync).toHaveBeenCalled();
    expect(scene.meshes.length).toBe(6);
  });
});

