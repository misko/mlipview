import { computeBonds } from '../../public/domain/bonding/computeBonds.js';

function vec(x, y, z) {
  return [x, y, z];
}

describe('computeBonds', () => {
  test('detects covalent pair in non-periodic space', () => {
    const result = computeBonds({
      elements: ['C', 'H'],
      positions: [vec(0, 0, 0), vec(1.08, 0, 0)],
    });
    expect(result.bonds).toHaveLength(1);
    const bond = result.bonds[0];
    expect(bond.i).toBe(0);
    expect(bond.j).toBe(1);
    expect(bond.crossing).toBe(false);
    expect(Array.isArray(bond.imageDelta)).toBe(true);
  });

  test('emits ghost info for periodic neighbours', () => {
    const result = computeBonds({
      elements: ['C', 'C'],
      positions: [vec(0.1, 0, 0), vec(2.9, 0, 0)],
      cell: {
        enabled: true,
        a: { x: 3, y: 0, z: 0 },
        b: { x: 0, y: 3, z: 0 },
        c: { x: 0, y: 0, z: 3 },
        originOffset: { x: 0, y: 0, z: 0 },
      },
    });
    expect(Array.isArray(result.ghostAtoms)).toBe(true);
    expect(result.ghostAtoms.length).toBeGreaterThan(0);
    expect(Array.isArray(result.ghostBondMeta)).toBe(true);
    expect(result.ghostBondMeta.length).toBeGreaterThan(0);
    expect(result.diagnostics?.periodic).toBeDefined();
  });

  test('spawns six ghost bonds for close pair under periodic boundary', () => {
    const positions = [vec(0.0, 0, 0), vec(1.4, 0, 0)];
    const cell = {
      enabled: true,
      a: { x: 10, y: 0, z: 0 },
      b: { x: 0, y: 10, z: 0 },
      c: { x: 0, y: 0, z: 10 },
      originOffset: { x: 0, y: 0, z: 0 },
    };

    const nonPeriodic = computeBonds({ elements: ['C', 'C'], positions });
    expect(nonPeriodic.bonds).toHaveLength(1);
    expect(nonPeriodic.ghostBondMeta).toHaveLength(0);

    const periodic = computeBonds({ elements: ['C', 'C'], positions, cell });
    expect(periodic.bonds).toHaveLength(1);
    expect(periodic.ghostBondMeta).toHaveLength(6);
    const offsets = periodic.ghostBondMeta
      .map((meta) => meta.shiftB.join(',') + '|' + meta.shiftA.join(','))
      .sort();
    expect(offsets).toEqual([
      '-1,0,0|-1,0,0',
      '0,-1,0|0,-1,0',
      '0,0,-1|0,0,-1',
      '0,0,1|0,0,1',
      '0,1,0|0,1,0',
      '1,0,0|1,0,0',
    ]);
  });

  test('produces ghost bonds only across the cell boundary for distant pair', () => {
    const positions = [vec(-4.5, 0, 0), vec(4.5, 0, 0)];
    const cell = {
      enabled: true,
      a: { x: 10, y: 0, z: 0 },
      b: { x: 0, y: 10, z: 0 },
      c: { x: 0, y: 0, z: 10 },
      originOffset: { x: -5, y: -5, z: -5 },
    };

    const nonPeriodic = computeBonds({ elements: ['C', 'C'], positions });
    expect(nonPeriodic.bonds).toHaveLength(0);
    expect(nonPeriodic.ghostBondMeta).toHaveLength(0);

    const periodic = computeBonds({ elements: ['C', 'C'], positions, cell });
    const crossingBonds = periodic.bonds.filter((b) => b.crossing);
    const primaryBonds = periodic.bonds.filter((b) => !b.crossing);
    expect(primaryBonds).toHaveLength(0);
    expect(crossingBonds).toHaveLength(1);
    expect(crossingBonds[0].opacity).toBe(0);
    expect(crossingBonds[0].imageDelta.some((v) => v !== 0)).toBe(true);
    expect(periodic.ghostBondMeta).toHaveLength(2);
    const offsets = periodic.ghostBondMeta
      .map((meta) => `${meta.shiftA.join(',')}|${meta.shiftB.join(',')}`)
      .sort();
    expect(offsets).toEqual(['-1,0,0|-1,0,0', '1,0,0|1,0,0']);
  });

  test('near-edge atoms generate two ghost bonds under periodic boundary', () => {
    const positions = [vec(9.5, 0, 0), vec(-9.5, 0, 0)];
    const cell = {
      enabled: true,
      a: { x: 10, y: 0, z: 0 },
      b: { x: 0, y: 10, z: 0 },
      c: { x: 0, y: 0, z: 10 },
      originOffset: { x: 0, y: 0, z: 0 },
    };

    const periodic = computeBonds({ elements: ['C', 'C'], positions, cell });
    const crossingBonds = periodic.bonds.filter((b) => b.crossing);
    const primaryBonds = periodic.bonds.filter((b) => !b.crossing);
    expect(primaryBonds).toHaveLength(0);
    expect(crossingBonds).toHaveLength(1);
    expect(crossingBonds[0].opacity).toBe(0);
    expect(crossingBonds[0].imageDelta.some((v) => v !== 0)).toBe(true);
    expect(periodic.ghostBondMeta).toHaveLength(2);
    const offsets = periodic.ghostBondMeta
      .map((meta) => `${meta.shiftA.join(',')}|${meta.shiftB.join(',')}`)
      .sort();
    expect(offsets).toEqual(['-1,0,0|-1,0,0', '1,0,0|1,0,0']);
  });
});

describe('computeBonds performance', () => {
  function createLattice(size) {
    const elements = [];
    const positions = [];
    for (let x = 0; x < size; x++) {
      for (let y = 0; y < size; y++) {
        for (let z = 0; z < size; z++) {
          elements.push((x + y + z) % 2 === 0 ? 'C' : 'H');
          positions.push(vec(x * 1.4, y * 1.4, z * 1.4));
        }
      }
    }
    return { elements, positions };
  }

  test('runs within expected time for 4³ lattice', () => {
    const { elements, positions } = createLattice(4); // 64 atoms
    const iterations = 3;
    const start = process.hrtime.bigint();
    for (let i = 0; i < iterations; i++) {
      computeBonds({
        elements,
        positions,
        cell: {
          enabled: true,
          a: { x: 10, y: 0, z: 0 },
          b: { x: 0, y: 10, z: 0 },
          c: { x: 0, y: 0, z: 10 },
          originOffset: { x: 0, y: 0, z: 0 },
        },
      });
    }
    const elapsedMs = Number(process.hrtime.bigint() - start) / 1e6;
    const average = elapsedMs / iterations;
    // Guard rail: keep average recompute under 1000ms on the CI host.
    expect(average).toBeLessThan(1000);
  });
});
