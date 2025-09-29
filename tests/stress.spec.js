import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createForceField } from '../public/physics/forcefield.js';
import { createDynamics } from '../public/physics/integrators.js';

describe('stress computation', () => {
  test('populates stress tensor after MD step', () => {
    const st = createMoleculeState({ elements:['C','C'], positions:[{x:0,y:0,z:0},{x:1.5,y:0,z:0}], bonds:[{i:0,j:1}] });
    createBondService(st);
    const ff = createForceField(st, {});
    const dyn = createDynamics(st, { timestep:0.2 });
    expect(st.dynamics.stress).toBeUndefined();
    dyn.stepMD({ forceFn: ff.computeForces, targetTemp:300 });
    expect(st.dynamics.stress).toBeDefined();
    expect(Number.isFinite(st.dynamics.stress.xx)).toBe(true);
  });
});
