import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';
import { createMoleculeView } from '../public/render/moleculeView.js';

// We rely on the jest BABYLON stub; no Engine/Scene constructor exists there, so create a minimal fake scene object.
function fakeScene() { return { onPointerObservable:{ add(){}, notify(){}, notifyObservers(){} } }; }

function buildSimple() {
  return { positions:[{x:0,y:0,z:0},{x:1.4,y:0,z:0}], elements:[6,6], bonds:[{i:0,j:1}] };
}

describe('highlight visibility', () => {
  test('selecting atom shows atom highlight mesh', () => {
  const mol = buildSimple();
  const state = createMoleculeState({ positions: mol.positions, elements: mol.elements, bonds: mol.bonds });
    const selection = createSelectionService(state);
    const scene = fakeScene();
    const view = createMoleculeView(scene, state);
    const { highlight } = view._internals;
    expect(highlight.atom.isVisible).toBe(false);
  selection.clickAtom(0);
    expect(highlight.atom.isVisible).toBe(true);
    expect(highlight.bond.isVisible).toBe(false);
  });
  test('selecting bond shows bond highlight mesh', () => {
  const mol = buildSimple();
  const state = createMoleculeState({ positions: mol.positions, elements: mol.elements, bonds: mol.bonds });
    const selection = createSelectionService(state);
    const scene = fakeScene();
    const view = createMoleculeView(scene, state);
    const { highlight } = view._internals;
  selection.clickBond({ i:0, j:1, key:'0-1', index:0 });
    expect(highlight.bond.isVisible).toBe(true);
    expect(highlight.atom.isVisible).toBe(false);
  });
});
