import { test, expect } from '@playwright/test';
import fs from 'fs';

const REF_PATH = 'public/reference/water_fairchem_bfgs_trace.json';
let REF = { energies: [] };
if (fs.existsSync(REF_PATH)) {
  try { REF = JSON.parse(fs.readFileSync(REF_PATH,'utf8')); } catch {}
}

async function serverReachable(page) {
  try {
    const resp = await page.request.post('http://127.0.0.1:8000/simple_calculate', { data: { atomic_numbers:[1], coordinates:[[0,0,0]] }});
    return resp.ok();
  } catch { return false; }
}

test('water FairChem BFGS 20-step first/last parity', async ({ page }) => {
  const up = await serverReachable(page);
  if (!up) {
    test.skip(true, 'FairChem server not reachable');
  }
  await page.goto('http://localhost:4000');
  await page.waitForFunction(()=>window.__MLIP_DEFAULT_LOADED === true);
  // Select water molecule
  await page.selectOption('.hud select:last-of-type','molecules/water.xyz');
  await page.waitForTimeout(200);
  // Switch provider to fairchem
  await page.selectOption('#forceProviderSel','fairchem');
  // Baseline already reset by provider switch; capture energies length
  await page.waitForFunction(()=> Array.isArray(window.__RELAX_TRACE) && window.__RELAX_TRACE.length === 1);
  // Perform 20 relax steps
  for (let i=0;i<20;i++) {
    await page.evaluate(()=>document.getElementById('btnRelax').click());
    await page.waitForFunction((expected)=> (window.__RELAX_TRACE||[]).length >= expected, i+2);
  }
  const energies = await page.evaluate(()=>window.__RELAX_TRACE);
  // Log first 21 energies (initial + 20 steps) for debugging/parity collection
  try {
    const first21 = energies.slice(0,21);
    // Structured prefix for easy grep
    console.log('E2E_FAIRCHEM_ENERGIES', JSON.stringify(first21));
  } catch(e) {}
  expect(energies.length).toBeGreaterThanOrEqual(21);
  const firstEnergy = energies[0];
  const lastEnergy = energies[energies.length-1];
  if (REF.energies && REF.energies.length>=2) {
    expect(Math.abs(firstEnergy - REF.energies[0])).toBeLessThan(5e-3);
    expect(Math.abs(lastEnergy - REF.energies[REF.energies.length-1])).toBeLessThan(5e-3);
  }
});
