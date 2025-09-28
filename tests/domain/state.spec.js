import { createMoleculeState, StateMutators, getDirtyFlags } from '../../public/domain/state/moleculeState.js';

// Basic tests for version bumps & dirty flags

describe('MoleculeState versions', () => {
  test('atoms and bonds bump versions independently', () => {
    const st = createMoleculeState();
    const prev = { ...st.versions };
    StateMutators.setAtoms(st, [{ id:0, element:'C', pos:{ x:0,y:0,z:0 } }]);
    expect(st.versions.atoms).toBe(prev.atoms + 1);
    expect(st.versions.bonds).toBe(prev.bonds);
    const dirty = getDirtyFlags(prev, st);
    expect(dirty.atoms).toBe(true);
    expect(dirty.bonds).toBeUndefined();
  });

  test('selection update emits event & bumps version', () => {
    const st = createMoleculeState();
    let called = 0;
    st.bus.on('selectionChanged', () => called++);
    StateMutators.setSelection(st, { kind:'atom', data:{ id:0 } });
    expect(st.versions.selection).toBe(1);
    expect(called).toBe(1);
  });
});
