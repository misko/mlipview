import { test, expect } from '@playwright/test';

test.setTimeout(30000);

async function loadWaterAndEnergy(page) {
  await page.waitForFunction(()=>window.__MLIP_DEFAULT_LOADED===true);
  // Use the dedicated molecule selector
  await page.selectOption('#moleculeSelect','molecules/water.xyz');
  await page.waitForTimeout(250);
  await page.evaluate(()=>{ window.__MLIPVIEW_SERVER='http://localhost:8000'; });
  // Recompute baseline energy synchronously to seed energySeries correctly
  await page.evaluate(async ()=>{ try { await window.viewerApi?.baselineEnergy?.(); } catch(e){} });
  // Poll computeForces until energy numeric
  await page.evaluate(async ()=>{
    for (let i=0;i<30;i++) {
      try { window.viewerApi?.ff?.computeForces?.(); } catch {}
      if (typeof window.viewerApi?.state?.dynamics?.energy === 'number' && isFinite(window.viewerApi.state.dynamics.energy)) break;
      await new Promise(r=>setTimeout(r,200));
    }
  });
  await page.waitForFunction(()=> typeof (window.viewerApi?.state?.dynamics?.energy) === 'number' && isFinite(window.viewerApi.state.dynamics.energy));
  return await page.evaluate(()=> window.viewerApi.state.dynamics.energy);
}

test('direct viewerApi.relaxStep() strict UMA energy baseline', async ({ page }) => {
  await page.addInitScript(()=>{ window.__FORCE_SYNC_REMOTE__ = true; window.__MLIPVIEW_SERVER='http://localhost:8000'; });
  await page.goto('/');
  const before = await loadWaterAndEnergy(page);
  // Geometry parity dump
  const geom = await page.evaluate(()=> window.__dumpCurrentAtoms && window.__dumpCurrentAtoms());
  expect(geom.atomic_numbers).toEqual([1,1,8]);
  // Strict energy expectation: large negative ~ -2079.
  expect(before).toBeLessThan(-1000);
  expect(before).toBeGreaterThan(-3000);
  await page.evaluate(async ()=>{ await window.viewerApi?.relaxStep?.(); });
  await page.waitForTimeout(300); // allow update
  const after = await page.evaluate(()=> window.viewerApi.state.dynamics.energy);
  console.log('PW_DIRECT_RELAX_BEFORE_AFTER', before, after);
  // Expect magnitude similar and improvement small (< 0.01 eV typical)
  expect(Math.abs(after - before)).toBeLessThan(0.02);
});
