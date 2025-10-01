import { test, expect } from '@playwright/test';

test.setTimeout(30000);

async function loadWater(page){
  await page.waitForFunction(()=>window.__MLIP_DEFAULT_LOADED===true);
  await page.selectOption('#moleculeSelect','molecules/water.xyz');
  await page.waitForTimeout(200);
  await page.evaluate(async ()=>{ try { window.__MLIPVIEW_SERVER='http://localhost:8000'; await window.viewerApi?.baselineEnergy?.(); window.viewerApi?.ff?.computeForces?.(); } catch{} });
  await page.waitForFunction(()=> typeof (window.viewerApi?.state?.dynamics?.energy) === 'number' && isFinite(window.viewerApi.state.dynamics.energy), null, { timeout: 15000 });
  await page.waitForFunction(()=> Array.isArray(window.__RELAX_TRACE) && window.__RELAX_TRACE.length===1);
}

test('UI relax step strict UMA energy baseline', async ({ page }) => {
  await page.addInitScript(()=>{ window.__FORCE_SYNC_REMOTE__=true; window.__MLIPVIEW_SERVER='http://localhost:8000'; });
  await page.goto('/');
  await loadWater(page);
  const before = await page.evaluate(()=> window.viewerApi?.state?.dynamics?.energy);
  const geom = await page.evaluate(()=> window.__dumpCurrentAtoms && window.__dumpCurrentAtoms());
  expect(geom.atomic_numbers).toEqual([1,1,8]);
  expect(before).toBeLessThan(-1000);
  expect(before).toBeGreaterThan(-3000);
  await page.click('#btnRelax');
  await page.waitForTimeout(300);
  const after = await page.evaluate(()=> window.viewerApi.state.dynamics.energy);
  console.log('PW_UI_RELAX_BEFORE_AFTER', before, after);
  expect(Math.abs(after - before)).toBeLessThan(0.02);
});
