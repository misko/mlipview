import { test, expect } from './fixtures.js';

test.setTimeout(30000);

async function loadWater(page) {
  const base = process.env.BASE_URL || 'http://127.0.0.1:5174';
  await page.goto(`${base}/?mol=molecules/water.xyz`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true);
  await page.evaluate(async () => {
    try {
      window.__MLIPVIEW_SERVER = 'http://localhost:8000';
      const ws = window.__fairchem_ws__ || window.__WS_API__;
      if (ws && typeof ws.ensureConnected === 'function') {
        await ws.ensureConnected();
        const st = window.viewerApi?.state;
        if (st) {
          const Z = (st.elements || []).map((e) => (typeof e === 'number' ? e : 0));
          const R = (st.positions || []).map((p) => [p.x, p.y, p.z]);
          ws.initSystem({ atomic_numbers: Z, positions: R });
        }
      }
      await window.viewerApi?.baselineEnergy?.();
      window.viewerApi?.ff?.computeForces?.();
    } catch {}
  });
  await page.waitForFunction(
    () =>
      typeof window.viewerApi?.state?.dynamics?.energy === 'number' &&
      isFinite(window.viewerApi.state.dynamics.energy),
    null,
    { timeout: 30000 }
  );
  await page.waitForFunction(
    () => Array.isArray(window.__RELAX_TRACE) && window.__RELAX_TRACE.length === 1
  );
}

test('UI relax step strict UMA energy baseline', async ({ page }) => {
  await page.addInitScript(() => {
    window.__FORCE_SYNC_REMOTE__ = true;
    window.__MLIPVIEW_SERVER = 'http://localhost:8000';
  });
  await page.goto('/');
  await loadWater(page);
  const before = await page.evaluate(() => window.viewerApi?.state?.dynamics?.energy);
  const geom = await page.evaluate(() => window.__dumpCurrentAtoms && window.__dumpCurrentAtoms());
  expect(geom.atomic_numbers).toEqual([1, 1, 8]);
  expect(before).toBeLessThan(-1000);
  expect(before).toBeGreaterThan(-3000);
  await page.click('#btnRelax');
  await page.waitForTimeout(300);
  const after = await page.evaluate(() => window.viewerApi.state.dynamics.energy);
  console.log('PW_UI_RELAX_BEFORE_AFTER', before, after);
  expect(Math.abs(after - before)).toBeLessThan(0.02);
});
