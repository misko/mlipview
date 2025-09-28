// SelectionService: wraps existing selection-model logic with a small facade.
import { createEmptySelection, applyBondClick, orientationToSide, nextBondCycleState } from '../../selection-model.js';

export function createSelectionService(state) {
  if (!state.selection || !('kind' in state.selection)) {
    state.selection = createEmptySelection();
  }
  return {
    clickBond(bondRef) { // bondRef: { i,j,key,index }
      const result = applyBondClick(state.selection, bondRef);
      state.versions.selection++; state.bus.emit('selectionChanged', { selection: state.selection, action: 'bondClick', result });
      return result;
    },
    clear() {
      state.selection = createEmptySelection();
      state.versions.selection++; state.bus.emit('selectionChanged', { selection: state.selection, action: 'clear' });
    },
    getSelection() { return state.selection; },
    helpers: { orientationToSide, nextBondCycleState }
  };
}
