import { createSimulationState } from '../public/physics/sim-model.js';
import { createCellOptimizer } from '../public/physics/cell-optimizer.js';

// Mock provider: stress diagonal proportional to (length - target)
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

describe('cell optimizer', () => {
  test('reduces cell lengths toward target isotropically', async () => {
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
    expect(Math.abs(after[0] - 10)).toBeLessThan(0.3);
    expect(res.converged).toBe(true);
  });

  test('no-op when provider returns no stress', async () => {
    const st = createSimulationState({
      Z: [6],
      positions: [[0, 0, 0]],
      cell: [
        [11, 0, 0],
        [0, 11, 0],
        [0, 0, 11],
      ],
    });
    const provider = {
      async compute() {
        return { energy: 0, forces: [], stress: null };
      },
    };
    const opt = createCellOptimizer(st, { forceProvider: provider, settings: { cellGamma: 0.5 } });
    const before = st.box.lengths[0];
    const r = await opt.step();
    expect(r.converged).toBe(true);
    expect(st.box.lengths[0]).toBe(before);
  });
});
