import { createEmptySelection } from '../public/selection-model.js';
import { createSelectionService } from '../public/domain/selectionService.js';

function makeMolState() {
  return {
    selection: createEmptySelection(),
    elements: Array.from({ length: 5 }, (_, i) => ({ Z: 6, index: i })),
    markSelectionChangedCalls: 0,
    markSelectionChanged() {
      this.markSelectionChangedCalls += 1;
    },
  };
}

describe('x-selection-null-guard', () => {
  test('clickAtom repopulates selection when state.selection is null', () => {
    const state = makeMolState();
    const svc = createSelectionService(state);
    state.selection = null;
    expect(() => svc.clickAtom(0)).not.toThrow();
    expect(state.selection).toBeTruthy();
    expect(state.markSelectionChangedCalls).toBeGreaterThan(0);
  });

  test('clickBond also tolerates null selection object', () => {
    const state = makeMolState();
    const svc = createSelectionService(state);
    state.selection = null;
    const bond = { i: 0, j: 1, key: '0-1', index: 0 };
    expect(() => svc.clickBond(bond)).not.toThrow();
    expect(state.selection).toBeTruthy();
  });
});
