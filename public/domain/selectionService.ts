import { createEmptySelection, applyBondClick } from '../selection-model.js';
import { computeBondRotationGroup } from './bondRotationUtils.js';
import { __count } from '../util/funcCount.ts';

export type BondRotationGroup = ReturnType<typeof computeBondRotationGroup>;

declare global {
  interface Window {
    __MLIPVIEW_DEBUG_SELECT?: boolean;
    __MLIP_LAST_SELECTION__?: unknown;
  }
}

export type BondSelectionData = {
  i: number;
  j: number;
  key?: unknown;
  index?: number;
  orientation?: number;
  rotationGroup?: BondRotationGroup;
};

export type BondSelection = { kind: 'bond'; data: BondSelectionData };
export type AtomSelection = { kind: 'atom'; data: { index: number; rotationGroup?: BondRotationGroup } };
export type NullSelection = { kind: null; data: null };
export type Selection = AtomSelection | BondSelection | NullSelection;

export interface SelectionState {
  elements: string[];
  selection?: Selection | null;
  markSelectionChanged: () => void;
}

type BondClickResult = ReturnType<typeof applyBondClick>;

function isSelection(value: unknown): value is Selection {
  if (!value || typeof value !== 'object') return false;
  const kind = (value as any).kind;
  return kind === 'atom' || kind === 'bond' || kind === null;
}

function ensureSelection(value: unknown): Selection {
  if (isSelection(value)) return value;
  return createEmptySelection() as Selection;
}

function makeEmptySelection(): Selection {
  return createEmptySelection() as Selection;
}

function isBondSelection(selection: Selection): selection is BondSelection {
  return selection.kind === 'bond';
}

function isAtomSelection(selection: Selection): selection is AtomSelection {
  return selection.kind === 'atom';
}

export function createSelectionService(molState: SelectionState) {
  __count('selectionService#createSelectionService');
  molState.selection = ensureSelection(molState.selection);

  function clearBondSpecific(): void {
    __count('selectionService#clearBondSpecific');
    molState.selection = ensureSelection(molState.selection);
    if (isBondSelection(molState.selection)) {
      molState.selection = makeEmptySelection();
    }
  }

  function clearAtomSpecific(): void {
    __count('selectionService#clearAtomSpecific');
    molState.selection = ensureSelection(molState.selection);
    if (isAtomSelection(molState.selection)) {
      molState.selection = makeEmptySelection();
    }
  }

  return {
    clickBond(ref: { i: number; j: number; key?: unknown; index?: number }): BondClickResult {
      __count('selectionService#clickBond');
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] clickBond', ref);
      } catch {}

      molState.selection = ensureSelection(molState.selection);

      if (isAtomSelection(molState.selection)) {
        molState.selection = makeEmptySelection();
      }

      const result = applyBondClick(molState.selection as any, ref);
      molState.selection = ensureSelection(molState.selection);

      if (isBondSelection(molState.selection)) {
        try {
          const { i, j, orientation } = molState.selection.data;
          molState.selection.data.rotationGroup = computeBondRotationGroup(
            molState as unknown as Parameters<typeof computeBondRotationGroup>[0],
            { i, j, orientation }
          );
        } catch {
          try {
            if (isBondSelection(molState.selection)) delete molState.selection.data.rotationGroup;
          } catch {}
        }
      }

      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] result bond', result, 'selection=', molState.selection);
      } catch {}
      molState.markSelectionChanged();
      return result;
    },

    clickAtom(index: number): 'ignored' | 'unchanged' | 'selected' {
      __count('selectionService#clickAtom');
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] clickAtom', index);
      } catch {}

      if (index == null || index < 0 || index >= molState.elements.length) return 'ignored';

      molState.selection = ensureSelection(molState.selection);
      if (isAtomSelection(molState.selection) && molState.selection.data.index === index) {
        return 'unchanged';
      }
      if (isBondSelection(molState.selection)) {
        molState.selection = makeEmptySelection();
      }

      molState.selection = { kind: 'atom', data: { index } };
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] result atom selection=', molState.selection);
      } catch {}
      molState.markSelectionChanged();
      try {
        if (typeof window !== 'undefined') {
          (window as any).__MLIP_LAST_SELECTION__ = { kind: 'atom', index };
        }
      } catch {}
      return 'selected';
    },

    clear(): void {
      __count('selectionService#clear');
      try {
        if (typeof window !== 'undefined' && window.__MLIPVIEW_DEBUG_SELECT)
          console.log('[selectService] clear');
      } catch {}
      molState.selection = makeEmptySelection();
      molState.markSelectionChanged();
    },

    get(): Selection {
      __count('selectionService#get');
      molState.selection = ensureSelection(molState.selection);
      return molState.selection;
    },

    _internal: { clearBondSpecific, clearAtomSpecific },
  };
}
