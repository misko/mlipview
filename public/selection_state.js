// public/selection_state.js
// Minimal shared selection coordinator to enforce mutual exclusivity
// between desktop atom selection and bond selection.

let kind = null; // 'atom' | 'bond' | null
let clearAtomFn = () => {};
let clearBondFn = () => {};
let atomDragging = false;
let _clearing = false; // reentrancy guard

export function registerAtomClear(fn) {
  if (typeof fn === 'function') clearAtomFn = fn;
}
export function registerBondClear(fn) {
  if (typeof fn === 'function') clearBondFn = fn;
}

export function selectAtom() {
  if (kind !== 'atom') {
    // Clear any bond selection when switching to atom
    if (!_clearing) {
      try { console.log('[select_state] selectAtom: clearing bond selection (prev kind=', kind, ')'); clearBondFn(); } catch {}
    }
  }
  kind = 'atom';
  console.log('[select_state] kind=atom');
}

export function selectBond() {
  if (kind !== 'bond') {
    // Clear any atom selection when switching to bond
    if (!_clearing) {
      try { console.log('[select_state] selectBond: clearing atom selection (prev kind=', kind, ')'); clearAtomFn(); } catch {}
    }
  }
  kind = 'bond';
  console.log('[select_state] kind=bond');
}

export function clearSelection() {
  if (_clearing) { console.log('[select_state] clearSelection ignored (reentrant)'); return; }
  _clearing = true;
  console.log('[select_state] clearSelection called (prev kind=', kind, ')');
  try { clearAtomFn(); } catch {}
  try { clearBondFn(); } catch {}
  kind = null;
  console.log('[select_state] kind=null');
  _clearing = false;
}

export function getKind() {
  return kind;
}

export function setAtomDragging(v) {
  atomDragging = !!v;
  console.log('[select_state] setAtomDragging', atomDragging);
}

export function isAtomDragging() {
  return atomDragging;
}
