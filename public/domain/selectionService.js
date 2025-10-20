import { createEmptySelection, applyBondClick } from '../selection-model.js';
import { __count } from '../util/funcCount.js';

export function createSelectionService(molState) {
  __count('selectionService#createSelectionService');
  if (!molState.selection || typeof molState.selection !== 'object') {
    molState.selection = createEmptySelection();
  }
  function clearBondSpecific() {
    __count('selectionService#clearBondSpecific');
    if (molState.selection.kind === 'bond') {
      molState.selection = createEmptySelection();
    }
  }
  function clearAtomSpecific() {
    __count('selectionService#clearAtomSpecific');
    if (molState.selection.kind === 'atom') {
      molState.selection = createEmptySelection();
    }
  }
  return {
    clickBond(ref) {
      __count('selectionService#clickBond');
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] clickBond', ref);
      } catch {}
      if (
        !molState.selection ||
        typeof molState.selection !== 'object' ||
        !('kind' in molState.selection)
      ) {
        molState.selection = createEmptySelection();
      }
      if (molState.selection.kind === 'atom') {
        molState.selection = createEmptySelection();
      }
      const result = applyBondClick(molState.selection, ref);
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] result bond', result, 'selection=', molState.selection);
      } catch {}
      molState.markSelectionChanged();
      return result;
    },
    clickAtom(index) {
      __count('selectionService#clickAtom');
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] clickAtom', index);
      } catch {}
      if (index == null || index < 0 || index >= molState.elements.length) return 'ignored';
      if (
        !molState.selection ||
        typeof molState.selection !== 'object' ||
        !('kind' in molState.selection)
      ) {
        molState.selection = createEmptySelection();
      }
      if (molState.selection.kind === 'atom' && molState.selection.data?.index === index)
        return 'unchanged';
      if (molState.selection.kind === 'bond') {
        molState.selection = createEmptySelection();
      }
      molState.selection = { kind: 'atom', data: { index } };
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] result atom selection=', molState.selection);
      } catch {}
      molState.markSelectionChanged();
      try {
        if (typeof window !== 'undefined') {
          window.__MLIP_LAST_SELECTION__ = { kind: 'atom', index };
        }
      } catch {}
      return 'selected';
    },
    clear() {
      __count('selectionService#clear');
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] clear');
      } catch {}
      molState.selection = createEmptySelection();
      molState.markSelectionChanged();
    },
    get() {
      __count('selectionService#get');
      return molState.selection;
    },
    _internal: { clearBondSpecific, clearAtomSpecific },
  };
}
