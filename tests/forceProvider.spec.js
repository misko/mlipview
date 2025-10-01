import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createLocalLJProvider } from '../public/physics/force-provider.js';

describe('force provider abstraction', () => {
  test('local provider returns forces and stress', async () => {
    const st = createMoleculeState({ elements:['C','C'], positions:[{x:0,y:0,z:0},{x:1.5,y:0,z:0}], bonds:[{i:0,j:1}] });
    createBondService(st);
    const provider = createLocalLJProvider(st, {});
    const res = await provider.compute({ wantStress:true });
    expect(Array.isArray(res.forces)).toBe(true);
    expect(res.forces.length).toBe(2);
    expect(res.forces[0].length).toBe(3);
    expect(res.stress).toBeTruthy();
    expect(Number.isFinite(res.stress.voigt[0])).toBe(true);
  });

  // FairChem HTTP provider tests removed after migration to server-side relax API only.
});
