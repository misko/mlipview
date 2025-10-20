import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

function makeLinearMolecule() {
  const state = createMoleculeState({
    elements: ['C', 'C', 'C', 'H', 'H'],
    positions: [
      { x: 0, y: 0, z: 0 }, // 0
      { x: 1, y: 0, z: 0 }, // 1
      { x: 2, y: 0, z: 0 }, // 2
      { x: 2.5, y: 1, z: 0 }, // 3 (attached to 2)
      { x: 2.5, y: -1, z: 0 }, // 4 (attached to 2)
    ],
    bonds: [
      { i: 0, j: 1 },
      { i: 1, j: 2 },
      { i: 2, j: 3 },
      { i: 2, j: 4 },
    ],
  });
  return state;
}

describe('manipulationService', () => {
  test('atom drag updates position and version', () => {
    const s = makeLinearMolecule();
    s.selection = { kind: 'atom', data: { index: 1 } };
    const manip = createManipulationService(s);
    const startVersion = s.versions.positions;
    manip.beginDrag(() => ({ x: 1, y: 0, z: 0 }));
    manip.updateDrag(() => ({ x: 2, y: 0, z: 0 })); // move to x=2
    expect(s.positions[1].x).toBeCloseTo(2, 5);
    expect(s.versions.positions).toBeGreaterThan(startVersion);
  });
  test('bond rotation rotates distal side', () => {
    const s = makeLinearMolecule();
    // select bond 1-2 orientation0 -> side 'j' (orientationToSide(0) => 'j') rotates the j side (atom index 2 subtree)
    s.selection = { kind: 'bond', data: { i: 1, j: 2, orientation: 0 } };
    const manip = createManipulationService(s);
    const preY3 = s.positions[3].y; // 1
    manip.rotateBond(Math.PI / 2); // rotate +90 around axis from j(2) anchored at i or anchor logic in service.
    const postY3 = s.positions[3].y;
    // y should change significantly (since rotation around axis should swing H)
    expect(Math.abs(postY3 - preY3)).toBeGreaterThan(0.1);
  });
});
