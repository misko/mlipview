import { createSimulationState } from '../public/physics/sim-model.js';
import { createOptimizer } from '../public/physics/optimizer.js';

// Mock force provider: harmonic well centered at origin with k=1
function createHarmonicProvider(k = 1) {
  return {
    capabilities: { supportsStress: false },
    async compute({ positions }) {
      let energy = 0;
      const forces = [];
      for (const p of positions) {
        const [x, y, z] = p;
        energy += 0.5 * k * (x * x + y * y + z * z);
        forces.push([-k * x, -k * y, -k * z]);
      }
      return { energy, forces };
    },
  };
}

describe('optimizer core', () => {
  test('steepest descent reduces distance to origin', async () => {
    const st = createSimulationState({
      Z: [6, 6],
      positions: [
        [2, 0, 0],
        [0, 2, 0],
      ],
    });
    const provider = createHarmonicProvider(1);
    const opt = createOptimizer(st, {
      forceProvider: provider,
      settings: { algorithm: 'sd', maxSteps: 200, gamma: 0.25, forceTol: 5e-4 },
    });
    const before = Math.hypot(st.pos[0], st.pos[1], st.pos[2]);
    const result = await opt.runUntil();
    const after = Math.hypot(st.pos[0], st.pos[1], st.pos[2]);
    expect(after).toBeLessThan(before);
    expect(result.converged).toBe(true);
  });
  test('conjugate gradient converges on harmonic', async () => {
    const st = createSimulationState({
      Z: [6, 6],
      positions: [
        [1.5, 0, 0],
        [0, 1.5, 0],
      ],
    });
    const provider = createHarmonicProvider(1);
    const opt = createOptimizer(st, {
      forceProvider: provider,
      settings: { algorithm: 'cg', maxSteps: 80, gamma: 0.3, forceTol: 1e-5, cgReset: 1000 },
    });
    const result = await opt.runUntil();
    expect(result.converged).toBe(true);
    // near origin
    expect(Math.abs(st.pos[0])).toBeLessThan(1e-2);
  });
});
