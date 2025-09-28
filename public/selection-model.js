// Shared selection model & helpers for desktop and VR
// A bond selection cycles: none -> orientation 0 (cyan) -> orientation 1 (purple) -> none

export function createEmptySelection() {
  return { kind: null, data: null }; // kind: 'bond' | 'atom' | null
}

export function isSameBond(sel, i, j) {
  return sel && sel.kind === 'bond' && sel.data && sel.data.i === i && sel.data.j === j;
}

export function nextBondCycleState(currentOrientation) {
  if (currentOrientation == null) return 0; // none -> 0
  if (currentOrientation === 0) return 1;   // 0 -> 1
  return null;                               // 1 -> none
}

export function applyBondClick(selection, bond) {
  // bond: { i, j, key, index }
  if (!selection) return null;
  if (isSameBond(selection, bond.i, bond.j)) {
    const next = nextBondCycleState(selection.data.orientation);
    if (next == null) {
      // Clear selection
      selection.kind = null;
      selection.data = null;
      return 'cleared';
    }
    selection.data.orientation = next;
    return next === 0 ? 'orientation0' : 'orientation1';
  }
  // New bond selection orientation 0
  selection.kind = 'bond';
  selection.data = { i: bond.i, j: bond.j, key: bond.key, index: bond.index, orientation: 0 };
  return 'orientation0';
}

export function clearSelection(selection) {
  if (!selection) return;
  selection.kind = null;
  selection.data = null;
}

export function orientationToSide(orientation) {
  // orientation 0 == original (i,j) => side 'j'
  // orientation 1 == reversed      => side 'i'
  return orientation === 1 ? 'i' : 'j';
}

export function bondOrientationColor(orientation) {
  if (orientation === 1) {
    return { diffuse: { r: 0.7, g: 0.4, b: 1.0 }, emissive: { r: 0.4, g: 0.15, b: 0.6 } };
  }
  // orientation 0 default cyan
  return { diffuse: { r: 0.1, g: 0.9, b: 1.0 }, emissive: { r: 0.05, g: 0.4, b: 0.5 } };
}

// Atom selection helper
export function selectAtom(selection, atom) { // atom: { idx, type }
  selection.kind = 'atom';
  selection.data = { idx: atom.idx, type: atom.type };
  return 'atom';
}
