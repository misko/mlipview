import { test, expect } from './fixtures.js';

test.setTimeout(30000);

async function loadWaterAndEnergy(page) {
  const base = process.env.BASE_URL || 'http://127.0.0.1:5174';
  await page.goto(`${base}/?mol=molecules/water.xyz`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true);
  await page.evaluate(() => {
    window.__MLIPVIEW_SERVER = 'http://localhost:8000';
  });
  // Recompute baseline energy synchronously to seed energySeries correctly
  await page.evaluate(async () => {
    try {
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
    } catch (e) {}
  });
  // Poll computeForces until energy numeric
  await page.evaluate(async () => {
    for (let i = 0; i < 30; i++) {
      try {
        window.viewerApi?.ff?.computeForces?.();
      } catch {}
      if (
        typeof window.viewerApi?.state?.dynamics?.energy === 'number' &&
        isFinite(window.viewerApi.state.dynamics.energy)
      )
        break;
      await new Promise((r) => setTimeout(r, 200));
    }
  });
  await page.waitForFunction(
    () =>
      typeof window.viewerApi?.state?.dynamics?.energy === 'number' &&
      isFinite(window.viewerApi.state.dynamics.energy),
    null,
    { timeout: 30000 }
  );
  return await page.evaluate(() => window.viewerApi.state.dynamics.energy);
}

test('direct viewerApi.relaxStep() strict UMA energy baseline', async ({ page }) => {
  await page.addInitScript(() => {
    window.__FORCE_SYNC_REMOTE__ = true;
    window.__MLIPVIEW_SERVER = 'http://localhost:8000';
  });
  await page.goto('/');
  const before = await loadWaterAndEnergy(page);
  // Geometry parity dump
  const geom = await page.evaluate(() => window.__dumpCurrentAtoms && window.__dumpCurrentAtoms());
  expect(geom.atomic_numbers).toEqual([1, 1, 8]);
  // Strict energy expectation: large negative ~ -2079.
  expect(before).toBeLessThan(-1000);
  expect(before).toBeGreaterThan(-3000);
  await page.evaluate(async () => {
    await window.viewerApi?.relaxStep?.();
  });
  await page.waitForTimeout(300); // allow update
  const after = await page.evaluate(() => window.viewerApi.state.dynamics.energy);
  console.log('PW_DIRECT_RELAX_BEFORE_AFTER', before, after);
  // Expect magnitude similar and improvement small (< 0.01 eV typical)
  expect(Math.abs(after - before)).toBeLessThan(0.02);
});
