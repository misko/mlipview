import fs from 'fs';
import http from 'http';
import { createBFGSStepper } from '../public/lj_bfgs.js';
import { fairchemCalculate } from '../public/fairchem_provider.js';

// Helper to check server availability quickly
function checkServer(url='http://127.0.0.1:8000') {
  return new Promise(resolve => {
    const req = http.request(url + '/simple_calculate', { method:'POST' }, res => { resolve(true); });
    req.on('error', ()=> resolve(false));
    req.setTimeout(500, ()=>{ try { req.destroy(); } catch{} resolve(false); });
    req.end(JSON.stringify({ atomic_numbers:[1], coordinates:[[0,0,0]] }));
  });
}

const REF_PATH = 'public/reference/water_fairchem_bfgs_trace.json';
let REF = { energies: [] };
if (fs.existsSync(REF_PATH)) {
  try { REF = JSON.parse(fs.readFileSync(REF_PATH,'utf8')); } catch {}
}

async function runStepperFairChem() {
  // Hardcode water geometry (same as water.xyz ordering H,H,O)
  const positions = [ [0.9575,0,0], [-0.2399872,0.92662721,0], [0,0,0] ];
  const elements = ['H','H','O'];
  const stepper = createBFGSStepper({ positions: positions.map(p=>p.slice()), fmax:0.05, maxStep:0.2, compute: async (p)=>{
    const { energy, forces } = await fairchemCalculate(p, elements);
    return { energy, forces };
  }});
  await stepper.getCurrent();
  for (let i=0;i<20;i++) await stepper.step();
  return stepper.history.map(h=>h.energy);
}

describe('fairchem water BFGS parity', () => {
  test('first & last energy vs reference (if available)', async () => {
    const up = await checkServer();
    if (!up) {
      console.warn('[fairchem parity] server not reachable, skipping test');
      return;
    }
    const energies = await runStepperFairChem();
    expect(energies.length).toBeGreaterThanOrEqual(2);
    if (REF.energies && REF.energies.length>=2) {
      const initDiff = Math.abs(energies[0] - REF.energies[0]);
      const finalDiff = Math.abs(energies[energies.length-1] - REF.energies[REF.energies.length-1]);
      expect(initDiff).toBeLessThan(5e-3);
      expect(finalDiff).toBeLessThan(5e-3);
    } else {
      console.warn('[fairchem parity] reference empty; skipping numeric asserts');
    }
  }, 30000);
});
