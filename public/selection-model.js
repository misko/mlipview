export function createEmptySelection() {
  return { kind: null, data: null };
}
export function isSameBond(sel, i, j) {
  return sel && sel.kind === 'bond' && sel.data && sel.data.i === i && sel.data.j === j;
}
export function nextBondCycleState(currentOrientation) {
  if (currentOrientation == null) return 0;
  if (currentOrientation === 0) return 1;
  return null;
}
export function applyBondClick(selection, bond) {
  if (!selection) return null;
  if (isSameBond(selection, bond.i, bond.j)) {
    const next = nextBondCycleState(selection.data.orientation);
    if (next == null) {
      selection.kind = null;
      selection.data = null;
      return 'cleared';
    }
    selection.data.orientation = next;
    return next === 0 ? 'orientation0' : 'orientation1';
  }
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
  return orientation === 1 ? 'i' : 'j';
}
export function bondOrientationColor(orientation) {
  if (orientation === 1) {
    return { diffuse: { r: 0.7, g: 0.4, b: 1.0 }, emissive: { r: 0.4, g: 0.15, b: 0.6 } };
  }
  return { diffuse: { r: 0.1, g: 0.9, b: 1.0 }, emissive: { r: 0.05, g: 0.4, b: 0.5 } };
}
export function selectAtom(selection, atom) {
  selection.kind = 'atom';
  selection.data = { idx: atom.idx, type: atom.type };
  return 'atom';
}
