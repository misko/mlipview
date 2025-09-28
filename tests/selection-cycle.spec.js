import { createEmptySelection, applyBondClick } from '../public/selection-model.js';

describe('selection model bond orientation cycle', () => {
  it('cycles orientation then clears', () => {
    const sel = createEmptySelection();
    // First click -> select bond orientation 0
  let r = applyBondClick(sel, { i: 1, j: 2, key: '1-2', index: 0 });
  expect(r).toBe('orientation0'); // first select => orientation0
    expect(sel.kind).toBe('bond');
    expect(sel.data.orientation).toBe(0);
    // Second click same bond -> orientation flips to 1
  r = applyBondClick(sel, { i: 1, j: 2, key: '1-2', index: 0 });
  expect(r).toBe('orientation1'); // second click => orientation1
    expect(sel.data.orientation).toBe(1);
    // Third click same bond -> clears
  r = applyBondClick(sel, { i: 1, j: 2, key: '1-2', index: 0 });
  expect(r).toBe('cleared'); // third click clears
    expect(sel.kind).toBeNull();
  });
});
