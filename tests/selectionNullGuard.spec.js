import { createEmptySelection } from '../public/selection-model.js';
import { createSelectionService } from '../public/domain/selectionService.js';

// Minimal molState stub mimicking required shape
function createMolState() {
  return {
    selection: createEmptySelection(),
    elements: new Array(5).fill(0).map((_, i) => ({ Z: 6, index: i })),
    markSelectionChangedCalls: 0,
    markSelectionChanged() {
      this.markSelectionChangedCalls++;
    },
  };
}

describe('selectionService null selection robustness', () => {
  test('clickAtom does not throw when molState.selection is null', () => {
    const molState = createMolState();
    const svc = createSelectionService(molState);
    // Simulate external code nulling selection (as seen in browser error)
    molState.selection = null; // <-- replicates bad state
    expect(() => svc.clickAtom(0)).not.toThrow();
    expect(molState.selection).toBeTruthy();
  });

  test('clickBond does not throw when molState.selection is null', () => {
    const molState = createMolState();
    const svc = createSelectionService(molState);
    molState.selection = null;
    const bondRef = { i: 0, j: 1, key: '0-1', index: 0 };
    expect(() => svc.clickBond(bondRef)).not.toThrow();
    expect(molState.selection).toBeTruthy();
  });
});
