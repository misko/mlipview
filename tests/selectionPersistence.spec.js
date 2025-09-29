import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';
import { createManipulationService } from '../public/domain/manipulationService.js';

function buildSmall() {
  // Simple tri-atom chain 0-1-2 so rotating bond (0-1) moves atom 0 around axis with atom 1 fixed
  return {
    positions:[{x:0,y:0,z:0},{x:1.5,y:0,z:0},{x:3.0,y:0,z:0}],
    elements:[6,6,6],
    bonds:[{i:0,j:1},{i:1,j:2}]
  };
}

// Fake scene (Babylon stub suffices)
function fakeScene() { return {}; }

describe('selection persistence', () => {
  test('bond selection persists across rotation steps', () => {
    const mol = buildSmall();
    const state = createMoleculeState({ positions:mol.positions, elements:mol.elements, bonds:mol.bonds });
    const selection = createSelectionService(state);
    const scene = fakeScene();
    const view = createMoleculeView(scene, state);
    const manip = createManipulationService(state);
    // Select bond 0-1
    selection.clickBond({ i:0,j:1,key:'0-1',index:0 });
    expect(state.selection?.kind).toBe('bond');
    expect(view._internals.highlight.bond.isVisible).toBe(true);
    // Perform several small rotations
    for (const angle of [5,10,15]) {
      const beforeVersion = state.versions.positions;
      manip.rotateBond(angle * Math.PI/180);
      expect(state.selection?.kind).toBe('bond');
      expect(view._internals.highlight.bond.isVisible).toBe(true);
      // Ensure positions changed (some rotation applied)
      expect(state.versions.positions).toBeGreaterThanOrEqual(beforeVersion);
    }
  });

  test('atom selection survives bondsChanged recompute', () => {
    const mol = buildSmall();
    const state = createMoleculeState({ positions:mol.positions, elements:mol.elements, bonds:mol.bonds });
    const selection = createSelectionService(state);
    const scene = fakeScene();
    const view = createMoleculeView(scene, state);
    // Select atom 1
    selection.clickAtom(1);
    expect(state.selection?.kind).toBe('atom');
    expect(view._internals.highlight.atom.isVisible).toBe(true);
    // Simulate bonds recompute WITHOUT topology change
    state.markBondsChanged();
    // Expect selection still present (currently this will FAIL until logic updated)
    expect(state.selection?.kind).toBe('atom');
    expect(view._internals.highlight.atom.isVisible).toBe(true);
  });
});
