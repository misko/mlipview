import { ljEnergyForces, bfgsOptimize } from '../public/lj_bfgs.js';

// Ground truth energies from ASE (provided): initial 2.802523, relaxed -2.983548 with fmax 0.05

function approx(a, b, tol) { return Math.abs(a - b) <= tol; }

describe('LJ + BFGS parity with ASE (water)', () => {
  test('initial and relaxed energies match reference within tolerance', async () => {
    const d = 0.9575;
    const t = Math.PI / 180 * 104.51;
    const positions = [
      [d, 0, 0],
      [d * Math.cos(t), d * Math.sin(t), 0],
      [0, 0, 0]
    ];
    const { energy: initialE } = ljEnergyForces(positions.map(p=>[...p]));
    // Tolerances: allow small numeric deviation due to float differences and missing neighbor list overhead
    expect(approx(initialE, 2.802523, 5e-3)).toBe(true);
    // Optimize
    const posCopy = positions.map(p=>[...p]);
    const result = await bfgsOptimize({ positions: posCopy, fmax: 0.05, maxSteps: 400, maxStep: 0.2, compute: async (p)=> ljEnergyForces(p) });
    expect(result.converged).toBe(true);
    const finalE = result.energy;
    expect(approx(finalE, -2.983548, 5e-3)).toBe(true);
  });
});
