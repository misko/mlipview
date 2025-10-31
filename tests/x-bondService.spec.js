import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.ts';

describe('x-bond service periodic', () => {
  test('detects bond across periodic boundary', () => {
    const st = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 2.9, y: 0, z: 0 },
      ],
      bonds: [],
    });
    st.cell = {
      a: { x: 3, y: 0, z: 0 },
      b: { x: 0, y: 3, z: 0 },
      c: { x: 0, y: 0, z: 3 },
      enabled: true,
      originOffset: { x: 0, y: 0, z: 0 },
    };

    const svc = createBondService(st);
    const bonds = svc.computePeriodicBonds();
    expect(bonds.some((b) => b.i === 0 && b.j === 1)).toBe(true);
  });
});
