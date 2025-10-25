import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createLocalLJProvider } from '../public/physics/force-provider.js';

describe('x-force-provider', () => {
  test('local LJ provider returns forces and stress', async () => {
    const state = createMoleculeState({
      elements: ['C', 'C'],
      positions: [
        { x: 0, y: 0, z: 0 },
        { x: 1.5, y: 0, z: 0 },
      ],
      bonds: [{ i: 0, j: 1 }],
    });
    createBondService(state);
    const provider = createLocalLJProvider(state, {});
    const res = await provider.compute({ wantStress: true });
    expect(Array.isArray(res.forces)).toBe(true);
    expect(res.forces.length).toBe(2);
    expect(res.stress).toBeTruthy();
  });
});

