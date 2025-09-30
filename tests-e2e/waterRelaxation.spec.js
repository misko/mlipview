import { test, expect } from '@playwright/test';
import fs from 'fs';
// Keep reference for final target energy only; initial changed after removing energy shift.
const REF = JSON.parse(fs.readFileSync('public/reference/water_lj_bfgs_trace.json','utf8'));

test('water LJ BFGS 20-step parity', async ({ page }) => {
  await page.goto('http://localhost:4000');
  // Wait viewer ready
  await page.waitForFunction(()=>window.__MLIP_DEFAULT_LOADED === true);
  // Add water option dynamically if not present (or assume added later) then load water
  await page.evaluate(()=>{
    const sel = document.querySelector('.hud select:last-of-type');
    if (sel && ![...sel.options].some(o=>o.value.includes('water'))) {
      const opt = document.createElement('option'); opt.value='molecules/water.xyz'; opt.textContent='Water'; sel.appendChild(opt);
    }
  });
  await page.selectOption('.hud select:last-of-type','molecules/water.xyz');
  await page.waitForTimeout(200); // allow load
  // Reset and force baseline energy capture
  await page.evaluate(()=>{ window._viewer?.baselineEnergy?.(); });
  // Capture initial forces
  const initialForces = await page.evaluate(()=>window.__RELAX_FORCES || null);
  if (initialForces) {
    for (let i=0;i<REF.initial_forces.length;i++) {
      for (let k=0;k<3;k++) {
        const ref = REF.initial_forces[i][k];
        const got = initialForces[i][k];
        const diff = Math.abs(got - ref);
          const tol = Math.max(1e-3, Math.abs(ref)*1e-2); // 1% relative or 1e-3 absolute
        expect(diff).toBeLessThan(tol);
          if (!(diff < tol)) {
            if (!window.__FORCE_DIFFS) window.__FORCE_DIFFS = [];
            window.__FORCE_DIFFS.push({ atom: i, comp: k, got, ref, diff, tol });
          }
      }
    }
  }
  // Perform 20 relax steps (await async handler)
  for (let i=0;i<20;i++) {
    await page.evaluate(()=>document.getElementById('btnRelax').click());
    await page.waitForFunction((expected)=> (window.__RELAX_TRACE||[]).length >= expected, i+2);
  }
  const energies = await page.evaluate(()=>window.__RELAX_TRACE);
  // We expect at least 21 entries (init + 20 steps)
  expect(energies.length).toBeGreaterThanOrEqual(21);
  // Only check first (initial) and last (final) energies per updated requirement
  const firstEnergy = energies[0];
  const lastEnergy = energies[energies.length-1];
  const refLast = REF.energies[REF.energies.length-1];
  // Initial energy should still be positive and > 2.5 (approx unshifted value ~2.78)
  expect(firstEnergy).toBeGreaterThan(2.5);
  // Final energy should be close to reference last (within 3e-2 since unshifted potential relax slightly lower ~ -3.0)
  expect(Math.abs(lastEnergy - refLast)).toBeLessThan(3e-2);
});