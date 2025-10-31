import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.ts';
import { createPickingService } from '../public/core/pickingService.js';

describe('x-atom drag camera', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
    BABYLON.PointerEventTypes = BABYLON.PointerEventTypes || { POINTERDOWN: 1 };
  });

  function makeScene({ hit = true } = {}) {
    return {
      pointerX: 50,
      pointerY: 50,
      pick: () => (hit ? { hit: true, pickedMesh: {}, thinInstanceIndex: 0 } : { hit: false }),
      onPointerObservable: {
        add(fn) {
          this._cb = fn;
        },
      },
      getEngine: () => ({
        getRenderingCanvas: () => ({ addEventListener() {} }),
      }),
    };
  }

  test('dragging atom leaves camera parameters unchanged', () => {
    const mol = createMoleculeState({
      elements: ['C', 'H'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    const selection = {
      get: () => mol.selection,
      clickAtom: (idx) => {
        mol.selection = { kind: 'atom', data: { index: idx } };
        mol.markSelectionChanged();
      },
      clickBond: () => {},
      clear: () => {
        mol.selection = { kind: null, data: null };
        mol.markSelectionChanged();
      },
    };
    const bondService = { recomputeAndStore: () => {} };
    const manipulation = createManipulationService(mol, { bondService });
    manipulation.molState = mol;

    const cameraState = {
      alpha: 1,
      beta: 2,
      radius: 10,
      target: { x: 0, y: 0, z: 0 },
    };
    const camera = {
      detachControl: () => {},
      attachControl: () => {},
      get alpha() {
        return cameraState.alpha;
      },
      get beta() {
        return cameraState.beta;
      },
      get radius() {
        return cameraState.radius;
      },
      get target() {
        return cameraState.target;
      },
    };

    const scene = makeScene({ hit: true });
    const view = {
      resolveAtomPick: () => ({ kind: 'atom', index: 0 }),
      resolveBondPick: () => null,
    };

    createPickingService(scene, view, selection, { manipulation, camera });
    selection.clickAtom(0);
    scene.onPointerObservable._cb?.({ type: BABYLON.PointerEventTypes.POINTERDOWN });

    const pre = {
      alpha: camera.alpha,
      beta: camera.beta,
      radius: camera.radius,
      target: { ...camera.target },
      position: { ...mol.positions[0] },
    };

    manipulation.updateDrag(() => ({ x: -0.5, y: 0, z: 0 }));

    expect(mol.positions[0].x).not.toBeCloseTo(pre.position.x, 6);
    expect(camera.alpha).toBe(pre.alpha);
    expect(camera.beta).toBe(pre.beta);
    expect(camera.radius).toBe(pre.radius);
    expect(camera.target.x).toBe(pre.target.x);
    expect(camera.target.y).toBe(pre.target.y);
    expect(camera.target.z).toBe(pre.target.z);
  });
});
