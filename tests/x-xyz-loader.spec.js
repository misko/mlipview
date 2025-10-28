import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('x-xyz-loader', () => {
  const xyz = `3
water
O 0 0 0
H 0.95 0 0
H -0.24 0.93 0`;

  test('parses elements and coordinates', () => {
    const parsed = parseXYZ(xyz);
    expect(parsed.elements).toEqual(['O', 'H', 'H']);
    expect(parsed.positions.length).toBe(3);
  });

  test('applies parsed data to molecule state', () => {
    const state = createMoleculeState();
    const parsed = parseXYZ(xyz);
    applyXYZToState(state, parsed);
    expect(state.elements[0]).toBe('O');
    expect(state.positions[2].x).toBeCloseTo(-0.24);
  });

  test('throws on malformed atom count line', () => {
    expect(() => parseXYZ('x\ncomment\n')).toThrow();
  });

  test('throws when atom lines are missing', () => {
    expect(() => parseXYZ('2\ncomment\nH 0 0 0')).toThrow();
  });
});

