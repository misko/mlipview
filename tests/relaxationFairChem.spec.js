import { createRelaxationRunner } from '../public/physics/relaxation-runner.js';
import { createSimulationState } from '../public/physics/sim-model.js';

// Mock FairChem-like provider returning diminishing forces + optional stress
function createMockFairChemProvider(){
  let scale = 1.0;
  return {
    async compute({ positions, wantStress }){
      const forces = positions.map(p=>{ const [x,y,z]=p; return [-0.1*scale*x,-0.1*scale*y,-0.1*scale*z]; });
      const energy = positions.reduce((e,p)=>{ const [x,y,z]=p; return e + 0.5*scale*(x*x+y*y+z*z); },0);
      scale *= 0.9; // decay forces
      let stress=null;
      if (wantStress) {
        const s = 0.05*scale;
        stress = { tensor:{ xx:s, yy:s, zz:s, xy:0,xz:0,yz:0 }, voigt:[s,s,s,0,0,0] };
      }
      return { energy, forces, stress };
    }
  };
}

describe('relaxation runner with mock FairChem', () => {
  test('converges with decaying forces and stress', async () => {
    const st = createSimulationState({ Z:[6,6], positions:[[1,0,0],[0,1,0]], cell:[[9,0,0],[0,9,0],[0,0,9]] });
    const provider = createMockFairChemProvider();
    const runner = createRelaxationRunner(st, { forceProvider: provider, settings:{ forceTol:1e-3, stressTol:1e-3, gamma:0.5, cellGamma:0.01, maxOuter:150, cellEvery:4 } });
    const res = await runner.run();
    expect(res.converged).toBe(true);
  });
});
