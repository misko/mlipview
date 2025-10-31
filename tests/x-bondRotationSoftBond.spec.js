import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.ts';

describe('x-bond rotation soft bond exclusion', () => {
  function makeChain() {
    const st = createMoleculeState({
      elements: ['C', 'C', 'C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1, y: 0, z: 0 },
        { x: 2, y: 0.5, z: 0 },
        { x: 3, y: 0.5, z: 0 },
      ],
      bonds: [
        { i: 0, j: 1 },
        { i: 1, j: 2 },
        { i: 2, j: 3 },
      ],
    });
    st.bonds = [
      { i: 0, j: 1, opacity: 1.0 },
      { i: 1, j: 2, opacity: 0.4 },
      { i: 2, j: 3, opacity: 1.0 },
    ];
    return st;
  }

  test('selected soft bond traverses across itself', () => {
    const st = makeChain();
    st.selection = { kind: 'bond', data: { i: 1, j: 2, orientation: 0 } };
    const manip = createManipulationService(st);
    manip.rotateBond(Math.PI / 10);
    const dbg = manip._debug.getLastRotation();
    expect(dbg).not.toBeNull();
    expect(dbg.sideAtoms.sort((a, b) => a - b)).toEqual([2, 3]);
  });

  test('soft non-selected bond blocks traversal', () => {
    const st = makeChain();
    st.selection = { kind: 'bond', data: { i: 0, j: 1, orientation: 0 } };
    const manip = createManipulationService(st);
    manip.rotateBond(Math.PI / 10);
    const dbg = manip._debug.getLastRotation();
    expect(dbg).not.toBeNull();
    expect(dbg.sideAtoms).toEqual([1]);
  });
});
