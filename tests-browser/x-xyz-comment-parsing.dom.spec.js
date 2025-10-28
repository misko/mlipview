/** @jest-environment jsdom */

import { parseXYZ } from '../public/util/xyzLoader.js';
import { getCellParameters } from '../public/util/pbc.js';

function xyzText({
  n = 2,
  comment = '',
  atoms = [
    ['H', 0, 0, 0],
    ['H', 0, 0, 0.74],
  ],
} = {}) {
  const lines = [String(n), comment];
  for (const a of atoms) {
    lines.push(`${a[0]} ${a[1]} ${a[2]} ${a[3]}`);
  }
  return lines.join('\n');
}

describe('XYZ comment parsing for temperature and cell', () => {
  test('extracts temperature and Lattice vectors', () => {
    const comment = 'foo temperature=350 Lattice=10,0,0;0,12,0;1,0,15';
    const txt = xyzText({ comment });
    const parsed = parseXYZ(txt);
    expect(parsed.temperature).toBe(350);
    expect(parsed.cell).toBeTruthy();
    const p = getCellParameters(parsed.cell);
    expect(p.a).toBeCloseTo(10, 6);
    expect(p.b).toBeCloseTo(12, 6);
    expect(p.c).toBeCloseTo(Math.hypot(1, 15), 6);
  });
  test('builds cell from abc+angles', () => {
    const comment = 'abc=10,12,15 alpha=90 beta=110 gamma=90 temp=298';
    const txt = xyzText({ comment });
    const parsed = parseXYZ(txt);
    expect(parsed.temperature).toBe(298);
    expect(parsed.cell).toBeTruthy();
    const p = getCellParameters(parsed.cell);
    expect(Math.round(p.beta)).toBe(110);
    expect(Math.round(p.alpha)).toBe(90);
    expect(Math.round(p.gamma)).toBe(90);
  });
});
