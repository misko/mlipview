import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createPickingService } from '../public/core/pickingService.js';

function fakeView() {
  return {
    resolveAtomPick(pick) {
      return pick?.hit ? { kind: 'atom', index: 0 } : null;
    },
    resolveBondPick() {
      return null;
    },
  };
}

describe('x-selection-clear', () => {
  test('empty pointer pick clears previous atom selection', () => {
    const state = createMoleculeState({
      elements: ['C'],
      positions: [{ x: 0, y: 0, z: 0 }],
    });
    const selection = createSelectionService(state);
    const scene = {
      pointerX: 0,
      pointerY: 0,
      pick: () => ({ hit: true }),
      onPointerObservable: { add() {} },
    };
    const view = fakeView();
    const picking = createPickingService(scene, view, selection);

    picking.pickAtPointer();
    selection.clickAtom(0);
    expect(selection.get().kind).toBe('atom');

    scene.pick = () => ({ hit: false });
    const res = picking.pickAtPointer();
    if (!res) selection.clear();
    expect(selection.get().kind).toBeNull();
  });
});
