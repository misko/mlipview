import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.js';

function makeState() {
  return createMoleculeState({ elements:['C','H','O'], positions:[{x:0,y:0,z:0},{x:1,y:0,z:0},{x:2,y:0,z:0}], bonds:[] });
}

describe('selectionService atom/bond mutual exclusivity', () => {
  test('select atom re-click is stable (no toggle)', () => {
    const st = makeState();
    const sel = createSelectionService(st);
    expect(sel.clickAtom(1)).toBe('selected');
    expect(sel.get().kind).toBe('atom');
    expect(sel.clickAtom(1)).toBe('unchanged');
    expect(sel.get().kind).toBe('atom');
  });
  test('bond clears atom selection first', () => {
    const st = makeState();
    st.bonds = [{ i:0,j:1 }];
    const sel = createSelectionService(st);
    sel.clickAtom(2);
    expect(sel.get().kind).toBe('atom');
    sel.clickBond({ i:0,j:1,key:'C-H',index:0 });
    expect(sel.get().kind).toBe('bond');
  });
  test('atom clears bond selection', () => {
    const st = makeState();
    st.bonds = [{ i:0,j:1 }];
    const sel = createSelectionService(st);
    sel.clickBond({ i:0,j:1,key:'C-H',index:0 });
    expect(sel.get().kind).toBe('bond');
    sel.clickAtom(0);
    expect(sel.get().kind).toBe('atom');
  });
});
