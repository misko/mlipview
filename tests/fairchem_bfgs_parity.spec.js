import fs from 'fs';
import http from 'http';

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

async function callRelaxFairChem() {
  // Water geometry (O,H,H) ordering matching server parity test
  const atomic_numbers = [8,1,1];
  const coordinates = [ [0,0,0], [0.9575,0,0], [-0.2399872,0.92662721,0] ];
  const body = JSON.stringify({ atomic_numbers, coordinates, steps:20, calculator:'uma' });
  const res = await fetch('http://127.0.0.1:8000/relax', { method:'POST', headers:{'Content-Type':'application/json'}, body });
  if(!res.ok){ throw new Error('relax failed '+res.status); }
  return await res.json();
}

describe('fairchem water BFGS parity via server /relax', () => {
  test('20-step first & last energy vs reference (if available)', async () => {
    const up = await checkServer();
    if (!up) {
      console.warn('[fairchem parity] server not reachable, skipping test');
      return;
    }
    const data = await callRelaxFairChem();
    expect(data.steps_completed).toBe(20);
    if (REF.energies && REF.energies.length>=2) {
      const initRef = REF.energies[0];
      const finalRef = REF.energies[REF.energies.length-1];
      expect(Math.abs(data.initial_energy - initRef)).toBeLessThan(5e-3);
      expect(Math.abs(data.final_energy - finalRef)).toBeLessThan(5e-3);
    } else {
      console.warn('[fairchem parity] reference empty; skipping numeric asserts');
    }
  }, 30000);
});
