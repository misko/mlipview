import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('XYZ loader', () => {
  const sample = `3\nwater\nO 0 0 0\nH 0.95 0 0\nH -0.24 0.93 0`;
  test('parses valid XYZ', () => {
    const parsed = parseXYZ(sample);
    expect(parsed.elements).toEqual(['O','H','H']);
    expect(parsed.positions.length).toBe(3);
  });
  test('applies to molecule state', () => {
    const st = createMoleculeState();
    const parsed = parseXYZ(sample);
    applyXYZToState(st, parsed);
    expect(st.elements[0]).toBe('O');
    expect(st.positions[2].x).toBeCloseTo(-0.24);
  });
  test('throws on malformed count', () => {
    expect(() => parseXYZ('x\ncomment\n')).toThrow();
  });
  test('throws on too few atom lines', () => {
    expect(() => parseXYZ('2\ncomment\nH 0 0 0')).toThrow();
  });
});
