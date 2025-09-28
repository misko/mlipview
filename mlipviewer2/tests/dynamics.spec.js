import { createMoleculeState } from '../public/domain/moleculeState.js';
import { createBondService } from '../public/domain/bondService.js';
import { createForceField } from '../public/physics/forcefield.js';
import { createDynamics } from '../public/physics/integrators.js';

describe('dynamics', () => {
  test('relax modifies positions', () => {
    const st = createMoleculeState({ elements:['C','C'], positions:[{x:0,y:0,z:0},{x:3,y:0,z:0}] });
    const bonds = createBondService(st); bonds.recomputeAndStore();
    const ff = createForceField(st, { r0:1.5 });
    const dyn = createDynamics(st, { timestep:0.1 });
    const before = st.positions[1].x;
    dyn.stepRelax({ forceFn: ff.computeForces, gamma: 0.01 });
    expect(st.positions[1].x).not.toBe(before);
  });
});
