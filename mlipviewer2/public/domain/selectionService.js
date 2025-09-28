import { createEmptySelection, applyBondClick } from './legacy/selection-model-adapter.js';

export function createSelectionService(molState) {
  if (!molState.selection || typeof molState.selection !== 'object') {
    molState.selection = createEmptySelection();
  }
  function clearBondSpecific() {
    if (molState.selection.kind === 'bond') {
      molState.selection = createEmptySelection();
    }
  }
  function clearAtomSpecific() {
    if (molState.selection.kind === 'atom') {
      molState.selection = createEmptySelection();
    }
  }
  return {
    clickBond(ref) {
      if (molState.selection.kind === 'atom') {
        molState.selection = createEmptySelection();
      }
      const result = applyBondClick(molState.selection, ref);
      molState.markSelectionChanged();
      return result;
    },
    clickAtom(index) {
      if (index == null || index < 0 || index >= molState.elements.length) return 'ignored';
      if (molState.selection.kind === 'atom' && molState.selection.data?.index === index) return 'unchanged';
      if (molState.selection.kind === 'bond') {
        molState.selection = createEmptySelection();
      }
      molState.selection = { kind: 'atom', data: { index } };
      molState.markSelectionChanged();
      return 'selected';
    },
    clear() {
      molState.selection = createEmptySelection();
      molState.markSelectionChanged();
    },
    get() { return molState.selection; },
    _internal: { clearBondSpecific, clearAtomSpecific }
  };
}
