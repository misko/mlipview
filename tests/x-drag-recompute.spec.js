import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.ts';
import { createBondService } from '../public/domain/bondService.ts';

function makeState() {
  return createMoleculeState({
    elements: ['C', 'H', 'H'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
      { x: -1, y: 0, z: 0 },
    ],
    bonds: [
      { i: 0, j: 1 },
      { i: 0, j: 2 },
    ],
  });
}

describe('x-drag-recompute', () => {
  test('dragging atom triggers bond recompute once', () => {
    const state = makeState();
    state.selection = { kind: 'atom', data: { index: 0 } };

    const bondService = createBondService(state);
    const spy = jest.spyOn(bondService, 'recomputeAndStore');
    const manip = createManipulationService(state, { bondService });

    const bondsVersion = state.versions.bonds;
    manip.beginDrag(() => ({ x: 0, y: 0, z: 0 }));
    manip.updateDrag(() => ({ x: 0.4, y: 0, z: 0 }));
    manip.endDrag();

    expect(state.versions.bonds).toBeGreaterThan(bondsVersion);
    expect(spy).toHaveBeenCalledTimes(1);
  });

  test('no movement leaves bond version unchanged', () => {
    const state = makeState();
    state.selection = { kind: 'atom', data: { index: 0 } };
    const bondService = createBondService(state);
    const spy = jest.spyOn(bondService, 'recomputeAndStore');
    const manip = createManipulationService(state, { bondService });

    const bondsVersion = state.versions.bonds;
    manip.beginDrag(() => ({ x: 0, y: 0, z: 0 }));
    manip.endDrag();

    expect(state.versions.bonds).toBe(bondsVersion);
    expect(spy).not.toHaveBeenCalled();
  });
});

