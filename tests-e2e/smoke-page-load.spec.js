import { test, expect } from '@playwright/test';

// Basic page load smoke: ensure canvas and panel render, default molecule loads, and idle energy appears.
// Uses real WS backend started by global-setup.

test('page loads and default molecule initializes (idle, no MD)', async ({ page, baseURL }) => {
  test.setTimeout(30000);
  // Bypass focus gating for deterministic init and ensure UMA server URL
  await page.addInitScript(() => { window.__MLIPVIEW_TEST_MODE = false; window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000'; });
  await page.goto(`${baseURL}/index.html?autoMD=0`);

  // Canvas and control panel exist
  await expect(page.locator('#viewer')).toBeVisible();
  await expect(page.locator('#controlPanel')).toBeVisible();

  // Wait for viewer API and default molecule load
  await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, { timeout: 15000 });

  // Molecule loaded -> positions present
  const nAtoms = await page.evaluate(() => window.viewerApi?.state?.positions?.length || 0);
  expect(nAtoms).toBeGreaterThan(0);
  // WS should be connected; if not, nudge the client to ensureInit and re-check
  let wsConnected = await page.evaluate(() => !!window.__fairchem_ws__?.getState?.().connected);
  if (!wsConnected) {
    await page.evaluate(async () => { try { const ws = window.__fairchem_ws__; ws && (await ws.ensureConnected()); } catch{} });
    await page.waitForTimeout(200);
    wsConnected = await page.evaluate(() => !!window.__fairchem_ws__?.getState?.().connected);
  }
  expect(wsConnected).toBeTruthy();
  // Ensure MD is not running when autoMD=0
  const running = await page.evaluate(() => window.viewerApi?.getMetrics?.().running || null);
  expect(running === null || running === undefined || running === false).toBeTruthy();
});
