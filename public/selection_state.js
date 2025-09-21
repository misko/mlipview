// public/selection_state.js
// Minimal shared selection coordinator to enforce mutual exclusivity
// between desktop atom selection and bond selection.

let kind = null; // 'atom' | 'bond' | null
let clearAtomFn = () => {};
let clearBondFn = () => {};
let atomDragging = false;

export function registerAtomClear(fn) {
  if (typeof fn === 'function') clearAtomFn = fn;
}
export function registerBondClear(fn) {
  if (typeof fn === 'function') clearBondFn = fn;
}

export function selectAtom() {
  if (kind !== 'atom') {
    // Clear any bond selection when switching to atom
    try { clearBondFn(); } catch {}
  }
  kind = 'atom';
}

export function selectBond() {
  if (kind !== 'bond') {
    // Clear any atom selection when switching to bond
    try { clearAtomFn(); } catch {}
  }
  kind = 'bond';
}

export function clearSelection() {
  try { clearAtomFn(); } catch {}
  try { clearBondFn(); } catch {}
  kind = null;
}

export function getKind() {
  return kind;
}

export function setAtomDragging(v) {
  atomDragging = !!v;
}

export function isAtomDragging() {
  return atomDragging;
}
