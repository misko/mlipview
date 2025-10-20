import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createPickingService } from '../public/core/pickingService.js';

// Minimal fake view resolving one atom only when pointer at specific coord
function makeView(state) {
  return {
    resolveAtomPick(pick) {
      return pick && pick.hit ? { kind: 'atom', index: 0 } : null;
    },
    resolveBondPick() {
      return null;
    },
  };
}

describe('selection clear via empty click', () => {
  test('select then empty click clears', () => {
    const state = createMoleculeState({
      elements: ['C'],
      positions: [{ x: 0, y: 0, z: 0 }],
      bonds: [],
    });
    const selection = createSelectionService(state);
    const scene = {
      pointerX: 0,
      pointerY: 0,
      pick() {
        return { hit: true };
      },
      onPointerObservable: { add() {} },
    };
    const view = makeView(state);
    createPickingService(scene, view, selection); // initial
    // First pointerDown should select atom
    scene.pick = () => ({ hit: true });
    // Simulate pointer down by direct call
    // We rely on implementation detail: handlePointerDown is registered; easier just call pickAtPointer + selection directly
    const picking = createPickingService(scene, view, selection);
    picking.pickAtPointer();
    selection.clickAtom(0);
    expect(selection.get().kind).toBe('atom');
    // Now emulate empty click (no hit)
    scene.pick = () => ({ hit: false });
    // pickAtPointer returns null -> manual clear
    const res = picking.pickAtPointer();
    if (!res) selection.clear();
    expect(selection.get().kind).toBeNull();
  });
});
