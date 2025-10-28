import { parseXYZ, applyXYZToState } from '../public/util/xyzLoader.js';
import { createMoleculeState } from '../public/domain/moleculeState.js';

describe('x-xyz-extended', () => {
  const xyzText = `3
Sample Lattice=1,0,0;0,2,0;0,0,3 temperature=300 note=TestCase
H 0 0 0
C 0.5 0.5 0.5
O 0.9 1.1 2.2`;

  test('parses lattice metadata and tags', () => {
    const parsed = parseXYZ(xyzText);
    expect(parsed.temperature).toBe(300);
    expect(parsed.tags.note).toBe('TestCase');
    expect(parsed.cell.a.x).toBe(1);
    expect(parsed.cell.b.y).toBe(2);
    expect(parsed.cell.c.z).toBe(3);
  });

  test('applies cell metadata to molecule state', () => {
    const parsed = parseXYZ(xyzText);
    const state = createMoleculeState();
    applyXYZToState(state, parsed);
    expect(state.cell.a.x).toBe(1);
    expect(state.cell.enabled).toBe(true);
  });
});

