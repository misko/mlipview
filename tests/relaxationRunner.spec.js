import { createSimulationState } from '../public/physics/sim-model.js';
import { createRelaxationRunner } from '../public/physics/relaxation-runner.js';

// Combined mock: harmonic atomic forces + elastic cell stress
function createCombinedProvider({ kAtom=1, targetLen=8, kCell=0.05 }){
  return {
    async compute({ positions, cell, wantStress }){
      let energy=0; const forces=[];
      for (const p of positions){ const [x,y,z]=p; energy += 0.5*kAtom*(x*x+y*y+z*z); forces.push([-kAtom*x, -kAtom*y, -kAtom*z]); }
      let stress=null;
      if (wantStress && cell){
        const aLen=Math.hypot(...cell[0]); const bLen=Math.hypot(...cell[1]); const cLen=Math.hypot(...cell[2]);
        const sxx=kCell*(aLen-targetLen); const syy=kCell*(bLen-targetLen); const szz=kCell*(cLen-targetLen);
        stress={ tensor:{ xx:sxx, yy:syy, zz:szz, xy:0,xz:0,yz:0 }, voigt:[sxx,syy,szz,0,0,0] };
      }
      return { energy, forces, stress };
    }
  };
}

describe('relaxation runner', () => {
  test('converges atomic + cell mock', async () => {
    const st = createSimulationState({ Z:[6,6], positions:[[2,0,0],[0,2,0]], cell:[[10,0,0],[0,10,0],[0,0,10]] });
    const provider = createCombinedProvider({ targetLen:8 });
    let updates=0; let last;
  const runner = createRelaxationRunner(st, { forceProvider: provider, settings:{ forceTol:5e-4, stressTol:5e-4, gamma:0.3, cellGamma:0.15, maxOuter:400, cellEvery:1, maxStrainStep:0.3, onUpdate:u=>{updates++; last=u;} } });
    const res = await runner.run();
    expect(res.converged).toBe(true);
    expect(updates).toBeGreaterThan(5);
    // Final cell length near target
  expect(Math.abs(st.box.lengths[0]-8)).toBeLessThan(0.6);
    // Atomic positions moved toward origin
    const r0 = Math.hypot(st.pos[0], st.pos[1], st.pos[2]);
    expect(r0).toBeLessThan(2);
  });
});
