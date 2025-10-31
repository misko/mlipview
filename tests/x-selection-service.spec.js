import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createSelectionService } from '../public/domain/selectionService.ts';

function makeState() {
  return createMoleculeState({
    elements: ['C', 'H', 'O'],
    positions: [
      { x: 0, y: 0, z: 0 },
      { x: 1, y: 0, z: 0 },
      { x: 2, y: 0, z: 0 },
    ],
    bonds: [{ i: 0, j: 1 }],
  });
}

describe('x-selection-service', () => {
  test('re-clicking same atom leaves selection unchanged', () => {
    const state = makeState();
    const svc = createSelectionService(state);
    expect(svc.clickAtom(1)).toBe('selected');
    expect(svc.clickAtom(1)).toBe('unchanged');
    expect(svc.get().kind).toBe('atom');
  });

  test('selecting bond after atom clears atom selection', () => {
    const state = makeState();
    const svc = createSelectionService(state);
    svc.clickAtom(2);
    expect(svc.get().kind).toBe('atom');
    svc.clickBond({ i: 0, j: 1, key: '0-1', index: 0 });
    expect(svc.get().kind).toBe('bond');
  });

  test('selecting atom clears prior bond selection', () => {
    const state = makeState();
    const svc = createSelectionService(state);
    svc.clickBond({ i: 0, j: 1, key: '0-1', index: 0 });
    expect(svc.get().kind).toBe('bond');
    svc.clickAtom(0);
    expect(svc.get().kind).toBe('atom');
  });
});
