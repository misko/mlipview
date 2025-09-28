// Jest tests for selection-model utilities
import { createEmptySelection, applyBondClick, orientationToSide, nextBondCycleState } from '../public/selection-model.js';

describe('selection-model bond cycle', () => {
  test('cycle none -> 0 -> 1 -> none', () => {
    const sel = createEmptySelection();
    expect(sel.kind).toBe(null);
    // first click
    let r = applyBondClick(sel, { i:0, j:1, key:'0-1', index:0 });
    expect(r).toBe('orientation0');
    expect(sel.data.orientation).toBe(0);
    // second click same bond -> orientation1
    r = applyBondClick(sel, { i:0, j:1, key:'0-1', index:0 });
    expect(r).toBe('orientation1');
    expect(sel.data.orientation).toBe(1);
    // third click clears
    r = applyBondClick(sel, { i:0, j:1, key:'0-1', index:0 });
    expect(r).toBe('cleared');
    expect(sel.kind).toBe(null);
  });

  test('select different bond resets cycle', () => {
    const sel = createEmptySelection();
    applyBondClick(sel, { i:0, j:1, key:'0-1', index:0 }); // orientation0
    applyBondClick(sel, { i:0, j:1, key:'0-1', index:0 }); // orientation1
    // new bond
    applyBondClick(sel, { i:2, j:3, key:'2-3', index:0 });
    expect(sel.kind).toBe('bond');
    expect(sel.data.i).toBe(2);
    expect(sel.data.orientation).toBe(0);
  });

  test('orientationToSide mapping', () => {
    expect(orientationToSide(0)).toBe('j');
    expect(orientationToSide(1)).toBe('i');
    expect(orientationToSide(undefined)).toBe('j');
  });

  test('nextBondCycleState helper', () => {
    expect(nextBondCycleState(undefined)).toBe(0);
    expect(nextBondCycleState(null)).toBe(0);
    expect(nextBondCycleState(0)).toBe(1);
    expect(nextBondCycleState(1)).toBe(null);
  });
});
