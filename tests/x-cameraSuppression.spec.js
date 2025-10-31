import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.ts';
import { createPickingService } from '../public/core/pickingService.js';

if (!global.BABYLON) global.BABYLON = {};
if (!BABYLON.PointerEventTypes) BABYLON.PointerEventTypes = { POINTERDOWN: 1 };
if (!BABYLON.Matrix) {
  BABYLON.Matrix = class {
    static Identity() {
      return new BABYLON.Matrix();
    }
  };
}

function makeSceneStub(pickResult) {
  const canvasListeners = {};
  const canvas = {
    addEventListener(type, fn) {
      canvasListeners[type] = fn;
    },
    getBoundingClientRect() {
      return { left: 0, top: 0 };
    },
  };
  return {
    pointerX: 0,
    pointerY: 0,
    pick() {
      return pickResult;
    },
    createPickingRay() {
      return {
        origin: {
          x: 0,
          y: 0,
          z: 0,
          add(v) {
            return { x: this.x + v.x, y: this.y + v.y, z: this.z + v.z };
          },
        },
        direction: {
          x: 0,
          y: 0,
          z: -1,
          scale() {
            return { x: 0, y: 0, z: -1 };
          },
        },
      };
    },
    onPointerObservable: {
      add(fn) {
        this._cb = fn;
      },
    },
    getEngine() {
      return {
        getRenderingCanvas() {
          return canvas;
        },
      };
    },
    __canvasListeners: canvasListeners,
  };
}

describe('x-camera suppression during drag', () => {
  test('detach/attach around drag lifecycle', () => {
    const st = createMoleculeState({
      elements: ['C'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    const selection = {
      clickAtom: (i) => {
        st.selection = { kind: 'atom', data: { index: i } };
        st.markSelectionChanged();
      },
      clickBond: () => {},
      clear: () => {},
      get: () => st.selection,
    };
    const manip = createManipulationService(st);

    const detachSpy = jest.fn();
    const attachSpy = jest.fn();
    const camera = { detachControl: detachSpy, attachControl: attachSpy };

    const scene = makeSceneStub({ hit: true, pickedMesh: {}, thinInstanceIndex: 0 });
    const view = {
      resolveAtomPick: () => ({ kind: 'atom', index: 0 }),
      resolveBondPick: () => null,
    };

    createPickingService(scene, view, selection, { manipulation: manip, camera });

    const canvasListeners = scene.__canvasListeners;
    canvasListeners.pointerdown?.({ clientX: 5, clientY: 5 });
    expect(detachSpy).toHaveBeenCalledTimes(1);

    canvasListeners.pointerup?.({ clientX: 5, clientY: 5 });
    expect(attachSpy).toHaveBeenCalledTimes(1);
  });
});
