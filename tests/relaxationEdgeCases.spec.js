import { createSimulationState } from '../public/physics/sim-model.js';
import { createRelaxationRunner } from '../public/physics/relaxation-runner.js';

// Provider returning zero forces & optional zero stress
function zeroProvider(){
  return { async compute({ positions, wantStress }) { return { energy:0, forces: positions.map(()=>[0,0,0]), stress: wantStress? { tensor:{xx:0,yy:0,zz:0,xy:0,xz:0,yz:0}, voigt:[0,0,0,0,0,0] }: null }; } };
}

// Provider with very small nearly-flat forces
function tinyForceProvider(scale=1e-6){
  return { async compute({ positions, wantStress }) { const forces = positions.map(p=>{ const [x,y,z]=p; return [-scale*x,-scale*y,-scale*z]; }); return { energy:0, forces, stress: wantStress? { tensor:{xx:scale,yy:scale,zz:scale,xy:0,xz:0,yz:0}, voigt:[scale,scale,scale,0,0,0] }: null }; } };
}

describe('relaxation edge cases', () => {
  test('already minimized (zero forces) converges immediately', async () => {
    const st = createSimulationState({ Z:[6], positions:[[0,0,0]] });
    const runner = createRelaxationRunner(st, { forceProvider: zeroProvider(), settings:{ maxOuter:10, forceTol:1e-5, stressTol:1e-5 } });
    const res = await runner.run();
    expect(res.converged).toBe(true);
    expect(res.steps).toBe(1);
  });
  test('nearly flat landscape still converges by force threshold', async () => {
    const st = createSimulationState({ Z:[6,6], positions:[[0.5,0,0],[0,0.5,0]] });
    const runner = createRelaxationRunner(st, { forceProvider: tinyForceProvider(), settings:{ maxOuter:50, forceTol:1e-5, stressTol:1e-4 } });
    const res = await runner.run();
    expect(res.converged).toBe(true);
  });
});
