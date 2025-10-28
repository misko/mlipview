import { createSimulationState } from '../public/physics/sim-model.js';
import { createCellOptimizer } from '../public/physics/cell-optimizer.js';

function createElasticCellProvider(targetLen = 10, k = 0.02) {
  return {
    async compute({ cell }) {
      const aLen = Math.hypot(cell[0][0], cell[0][1], cell[0][2]);
      const bLen = Math.hypot(cell[1][0], cell[1][1], cell[1][2]);
      const cLen = Math.hypot(cell[2][0], cell[2][1], cell[2][2]);
      const sxx = k * (aLen - targetLen);
      const syy = k * (bLen - targetLen);
      const szz = k * (cLen - targetLen);
      return {
        energy: 0,
        forces: [],
        stress: {
          tensor: { xx: sxx, yy: syy, zz: szz, xy: 0, xz: 0, yz: 0 },
          voigt: [sxx, syy, szz, 0, 0, 0],
        },
      };
    },
  };
}

describe('x-cell-optimizer (xfail)', () => {
  test.failing('reduces lattice lengths toward target (unsupported)', async () => {
    const st = createSimulationState({
      Z: [6],
      positions: [[0, 0, 0]],
      cell: [
        [12, 0, 0],
        [0, 12, 0],
        [0, 0, 12],
      ],
    });
    const provider = createElasticCellProvider(10, 0.05);
    const opt = createCellOptimizer(st, {
      forceProvider: provider,
      settings: {
        isotropic: true,
        cellGamma: 0.5,
        stressTol: 1e-3,
        maxSteps: 50,
        maxStrainStep: 0.2,
      },
    });
    const before = st.box.lengths.slice();
    const res = await opt.runUntil();
    const after = st.box.lengths;
    expect(after[0]).toBeLessThan(before[0]);
    expect(after[1]).toBeLessThan(before[1]);
    expect(after[2]).toBeLessThan(before[2]);
    // TODO: tighten this to a realistic tolerance once variable cell support returns.
    expect(Math.abs(after[0] - 10)).toBeLessThan(1e-3);
    expect(Math.abs(after[1] - 10)).toBeLessThan(1e-3);
    expect(Math.abs(after[2] - 10)).toBeLessThan(1e-3);
    expect(res.converged).toBe(true);
  });
});
