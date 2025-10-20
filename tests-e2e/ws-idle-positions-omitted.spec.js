import { test, expect } from '@playwright/test';

// Verify that during idle (no simulation running), idle compute frames omit positions
// and include energy (and likely forces). We hook the client WS decoded frames.

async function collectFrames(page, { timeout = 8000 } = {}) {
  return await page.evaluate(
    ({ timeout }) => {
      return new Promise((resolve) => {
        const out = [];
        const ws = window.__fairchem_ws__;
        if (!ws) return resolve(out);
        const off = ws.onResult((r) => {
          try {
            out.push(r);
          } catch {}
        });
        setTimeout(() => {
          try {
            off && off();
          } catch {}
          resolve(out);
        }, timeout);
      });
    },
    { timeout }
  );
}

test('idle frames omit positions and carry energy', async ({ page, baseURL }) => {
  test.setTimeout(45000);
  await page.addInitScript(() => {
    window.__MLIPVIEW_TEST_MODE = false;
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  });
  await page.goto(`${baseURL}/index.html?autoMD=0`);

  await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
    timeout: 45000,
  });
  // Trigger an idle compute explicitly
  await page.evaluate(() => window.viewerApi.requestSimpleCalculateNow());

  const frames = await collectFrames(page, { timeout: 3000 });
  const idle = frames.find((f) => typeof f?.energy === 'number');
  expect(idle).toBeTruthy();
  // Idle frames should not contain positions per protocol
  if (idle) {
    const hasPositions = Array.isArray(idle.positions) && idle.positions.length > 0;
    expect(hasPositions).toBeFalsy();
  }
});
