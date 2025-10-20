import {
  computeOrthoCellFromPositions,
  wrapPositionsInPlace,
  cellToMatrixArray,
} from '../public/util/pbc.js';

describe('PBC utilities', () => {
  test('computeOrthoCellFromPositions gives 1Å padding each side', () => {
    const pos = [
      { x: 0, y: 0, z: 0 },
      { x: 2, y: 3, z: 4 },
    ];
    const cell = computeOrthoCellFromPositions(pos, 1.0);
    expect(cell.enabled).toBe(true);
    expect(cell.synthetic).toBe(true);
    // Extents should be (max-min)+2*pad
    expect(cell.a.x).toBeCloseTo(2 - 0 + 2 * 1, 6);
    expect(cell.b.y).toBeCloseTo(3 - 0 + 2 * 1, 6);
    expect(cell.c.z).toBeCloseTo(4 - 0 + 2 * 1, 6);
    // Origin offset should be min - pad
    expect(cell.originOffset.x).toBeCloseTo(-1, 6);
    expect(cell.originOffset.y).toBeCloseTo(-1, 6);
    expect(cell.originOffset.z).toBeCloseTo(-1, 6);
  });

  test('wrapPositionsInPlace wraps into primary cell (orthorhombic)', () => {
    const pos = [
      { x: -2, y: 0.5, z: 0.5 },
      { x: 3.1, y: -0.2, z: 4.9 },
    ];
    const cell = {
      a: { x: 3, y: 0, z: 0 },
      b: { x: 0, y: 2, z: 0 },
      c: { x: 0, y: 0, z: 5 },
      originOffset: { x: 0, y: 0, z: 0 },
      enabled: true,
    };
    const ok = wrapPositionsInPlace(pos, cell);
    expect(ok).toBe(true);
    // After wrap: 0<=x<3, 0<=y<2, 0<=z<5
    for (const p of pos) {
      expect(p.x).toBeGreaterThanOrEqual(0 - 1e-9);
      expect(p.x).toBeLessThan(3 + 1e-9);
      expect(p.y).toBeGreaterThanOrEqual(0 - 1e-9);
      expect(p.y).toBeLessThan(2 + 1e-9);
      expect(p.z).toBeGreaterThanOrEqual(0 - 1e-9);
      expect(p.z).toBeLessThan(5 + 1e-9);
    }
  });

  test('cellToMatrixArray returns 3x3', () => {
    const cell = {
      a: { x: 1, y: 0.1, z: 0 },
      b: { x: 0, y: 2, z: 0.2 },
      c: { x: 0.3, y: 0, z: 3 },
    };
    const M = cellToMatrixArray(cell);
    expect(M.length).toBe(3);
    expect(M[0].length).toBe(3);
  });
});
