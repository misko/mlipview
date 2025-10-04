import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

/*
 * Tests rotation behavior with soft bonds (low opacity) excluded from traversal except the selected bond.
 */

describe('bond rotation exclusion of soft bonds', () => {
  function makeChain() {
    const st = createMoleculeState({
      elements:['C','C','C','C'],
      positions:[
        {x:0,y:0,z:0},   // 0
        {x:1,y:0,z:0},   // 1
        {x:2,y:0.5,z:0}, // 2 (offset)
        {x:3,y:0.5,z:0}  // 3
      ],
      bonds:[ {i:0,j:1}, {i:1,j:2}, {i:2,j:3} ]
    });
    st.bonds = [
      { i:0,j:1, opacity:1.0 },
      { i:1,j:2, opacity:0.4 }, // soft
      { i:2,j:3, opacity:1.0 }
    ];
    return st;
  }

  test('selected soft bond still allows traversal beyond moving root', () => {
    const st = makeChain();
    st.selection = { kind:'bond', data:{ i:1, j:2, orientation:0 } }; // bond 1-2 selected
    const manip = createManipulationService(st);
    manip.rotateBond(Math.PI/10);
    const dbg = manip._debug.getLastRotation();
    expect(dbg).not.toBeNull();
    // Expect sideAtoms includes movingRoot (2) and continues to 3 via strong bond 2-3
    expect(dbg.sideAtoms.sort((a,b)=>a-b)).toEqual([2,3]);
  });

  test('soft non-selected bond blocks traversal', () => {
    const st = makeChain();
    st.selection = { kind:'bond', data:{ i:0, j:1, orientation:0 } }; // bond 0-1 selected
    const manip = createManipulationService(st);
    manip.rotateBond(Math.PI/10);
    const dbg = manip._debug.getLastRotation();
    expect(dbg).not.toBeNull();
    // Expect sideAtoms only includes atom 1 (cannot cross soft 1-2)
    expect(dbg.sideAtoms).toEqual([1]);
  });
});
