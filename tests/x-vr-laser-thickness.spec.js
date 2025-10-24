import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';

describe('x-vr laser thickness regression', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
  });

  test('laser diameter grows on hover and relaxes afterwards', async () => {
    const pickSequence = [];
    for (let i = 0; i < 4; i++) pickSequence.push(null);
    for (let i = 0; i < 6; i++) pickSequence.push({ kind: 'atom', dist: 4 });
    for (let i = 0; i < 6; i++) pickSequence.push(null);
    let pickIndex = 0;

    const controllerState = new Map();
    const controllerNode = {
      getAbsolutePosition: () => new BABYLON.Vector3(0, 0, 0),
      getDirection: () => BABYLON.Vector3.Forward(),
    };
    const controller = {
      uniqueId: 'laser-regression',
      pointer: controllerNode,
      motionController: { getComponent: () => ({ pressed: false }) },
      inputSource: { gamepad: { buttons: [{ pressed: false }], axes: [] } },
    };
    const controllers = [controller];

    const scene = { meshes: [] };

    const { createVRLaserManager } = await import('../public/vr/vr-laser.js');
    const laserMgr = createVRLaserManager({
      scene,
      xrHelper: { input: { controllers } },
      controllerState,
      getControllerNode: () => controllerNode,
      pickAtomOrBondWithRay: () => pickSequence[pickIndex++] || null,
    });

    const step = () => {
      laserMgr.updateFrame(controllers);
      return laserMgr.getInfo()[0]?.diameter;
    };

    let baseDiam;
    for (let i = 0; i < 4; i++) baseDiam = step();
    expect(baseDiam).toBeDefined();

    let hoverDiam;
    for (let i = 0; i < 6; i++) hoverDiam = step();
    expect(hoverDiam).toBeGreaterThan(baseDiam);

    let postDiam;
    for (let i = 0; i < 6; i++) postDiam = step();
    expect(postDiam).toBeLessThan(hoverDiam);
  });
});
