import { jest } from '@jest/globals';
import { installBabylonStub } from './helpers/installBabylonStub.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.ts';
import { createPickingService } from '../public/core/pickingService.js';

describe('x-atom drag screen plane', () => {
  beforeEach(() => {
    jest.resetModules();
    installBabylonStub();
    BABYLON.PointerEventTypes = BABYLON.PointerEventTypes || { POINTERDOWN: 1 };
  });

  function setup() {
    const state = createMoleculeState({
      elements: ['C'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    const selection = {
      get: () => state.selection,
      clickAtom: (idx) => {
        state.selection = { kind: 'atom', data: { index: idx } };
        state.markSelectionChanged();
      },
      clickBond: () => {},
      clear: () => {
        state.selection = { kind: null, data: null };
        state.markSelectionChanged();
      },
    };
    const manipulation = createManipulationService(state, {});
    manipulation.molState = state;
    const camera = {
      position: new BABYLON.Vector3(0, 0, -10),
      detachControl: () => {},
      attachControl: () => {},
      get alpha() {
        return 0;
      },
      get beta() {
        return 0;
      },
      get radius() {
        return 10;
      },
      get target() {
        return new BABYLON.Vector3(0, 0, 0);
      },
    };

    const viewport = { cx: 100, cy: 100 };
    const makeRay = (x, y) => {
      const nx = (x - viewport.cx) / viewport.cx;
      const ny = (y - viewport.cy) / viewport.cy;
      const dir = new BABYLON.Vector3(nx, -ny, 1).normalize();
      return { origin: camera.position.clone(), direction: dir };
    };

    const scene = {
      pointerX: 100,
      pointerY: 100,
      pick: () => ({ hit: true, pickedMesh: {}, thinInstanceIndex: 0 }),
      onPointerObservable: {
        add(fn) {
          this._cb = fn;
        },
      },
      getEngine: () => ({
        getRenderingCanvas: () => ({ addEventListener() {}, removeEventListener() {} }),
      }),
      createPickingRay: (x, y) => makeRay(x, y),
    };

    const view = {
      resolveAtomPick: () => ({ kind: 'atom', index: 0 }),
      resolveBondPick: () => null,
    };

    createPickingService(scene, view, selection, { manipulation, camera });
    selection.clickAtom(0);
    scene.onPointerObservable._cb?.({ type: BABYLON.PointerEventTypes.POINTERDOWN });

    return { state, manipulation, scene };
  }

  test('drag plane normal aligns with camera vector', () => {
    const { manipulation } = setup();
    const debug = manipulation._debug.getDragState();
    expect(debug).toBeTruthy();
    const planeNormal = { ...debug.planeNormal };
    expect(Math.hypot(planeNormal.x, planeNormal.y, planeNormal.z)).toBeCloseTo(1, 6);
    expect(planeNormal.z).toBeGreaterThan(0);
  });

  test('pointer moves adjust atom within camera-aligned plane', () => {
    const { state, manipulation, scene } = setup();

    const planePoint = { x: 0, y: 0, z: 0 };
    const planeNormal = { x: 0, y: 0, z: 1 };

    const intersect = () => {
      const ray = scene.createPickingRay(scene.pointerX, scene.pointerY);
      const denom =
        planeNormal.x * ray.direction.x +
        planeNormal.y * ray.direction.y +
        planeNormal.z * ray.direction.z;
      const vx = planePoint.x - ray.origin.x;
      const vy = planePoint.y - ray.origin.y;
      const vz = planePoint.z - ray.origin.z;
      const t =
        (planeNormal.x * vx + planeNormal.y * vy + planeNormal.z * vz) /
        denom;
      return {
        x: ray.origin.x + ray.direction.x * t,
        y: ray.origin.y + ray.direction.y * t,
        z: ray.origin.z + ray.direction.z * t,
      };
    };

    scene.pointerX = 140;
    manipulation.updateDrag(() => intersect());

    expect(state.positions[0].x).toBeGreaterThan(0);
    expect(state.positions[0].z).toBeCloseTo(0, 6);

    scene.pointerY = 60;
    manipulation.updateDrag(() => intersect());

    expect(state.positions[0].y).toBeGreaterThan(0);
    expect(state.positions[0].z).toBeCloseTo(0, 6);
  });
});
