import { createMoleculeState } from '../../public/domain/state/moleculeState.js';
import { createSelectionService } from '../../public/domain/services/selectionService.js';

// Provide minimal selection-model dependency already tested elsewhere; here we integration-test facade.

describe('SelectionService facade', () => {
  test('bond click cycles orientations', () => {
    const st = createMoleculeState({ selection: { kind:null, data:{} } });
    const svc = createSelectionService(st);
    const r1 = svc.clickBond({ i:0, j:1, key:'0-1', index:0 });
    const r2 = svc.clickBond({ i:0, j:1, key:'0-1', index:0 });
    const r3 = svc.clickBond({ i:0, j:1, key:'0-1', index:0 });
    expect([r1,r2,r3]).toEqual(['orientation0','orientation1','cleared']);
  });

  test('clear resets selection', () => {
    const st = createMoleculeState();
    const svc = createSelectionService(st);
    svc.clickBond({ i:0, j:1, key:'0-1', index:0 });
    svc.clear();
    expect(svc.getSelection().kind).toBe(null);
  });
});
