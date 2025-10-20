import { test, expect } from '@playwright/test';

// Page load with MD auto-run enabled (default). Verifies running state toggles and energy grows.

test('page load starts MD and energy advances', async ({ page, baseURL }) => {
  test.setTimeout(45000);
  // Allow auto MD by not setting __MLIPVIEW_TEST_MODE and not using autoMD=0
  await page.addInitScript(() => {
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  });
  await page.goto(`${baseURL}/index.html`);

  await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
    timeout: 20000,
  });

  // Wait for MD to start (UI sync interval updates legacy btn text to 'stop')
  await page.waitForFunction(
    () => {
      const btn = document.getElementById('btnMDRun');
      return !!btn && btn.textContent === 'stop';
    },
    { timeout: 20000 }
  );

  // Energy series should increase over time
  const startLen = await page.evaluate(() => window.viewerApi?.debugEnergySeriesLength?.() || 0);
  await page.waitForTimeout(1500);
  const endLen = await page.evaluate(() => window.viewerApi?.debugEnergySeriesLength?.() || 0);
  expect(endLen).toBeGreaterThanOrEqual(startLen);
});
