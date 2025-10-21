// Purpose: End-to-end: an idle USER_INTERACTION echoes the frontend-provided
// user_interaction_count (UIC) in the server result, validating coalescing.
import { test, expect } from '@playwright/test';

// Verify user_interaction_count (UIC) is echoed by the server in idle computes
// when the client sends USER_INTERACTION with setCounters.

async function onceEnergy(page, { timeout = 5000 } = {}) {
  return await page.evaluate(
    ({ timeout }) => {
      const ws = window.__fairchem_ws__;
      return new Promise((resolve, reject) => {
        if (!ws) return resolve(null);
        const off = ws.onResult((r) => {
          if (r && typeof r.energy === 'number') {
            try {
              off && off();
            } catch { }
            resolve(r);
          }
        });
        setTimeout(() => {
          try {
            off && off();
          } catch { }
          reject(new Error('timeout'));
        }, timeout);
      });
    },
    { timeout }
  );
}

test('idle compute echoes user_interaction_count', async ({ page, baseURL }) => {
  test.setTimeout(45000);
  await page.addInitScript(() => {
    window.__MLIPVIEW_TEST_MODE = false;
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  });
  await page.goto(`${baseURL}/index.html?autoMD=0&wsDebug=1&debug=1`);
  await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
    timeout: 45000,
  });

  // Attach listener BEFORE sending USER_INTERACTION to avoid race and wait for matching UIC
  const echoed = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    const api = window.viewerApi;
    await ws.ensureConnected();
    // Choose a distinctive UIC value to match
    const cur = 987654321;
    return await new Promise((resolve) => {
      const off = ws.onResult((r) => {
        if (r && typeof r.energy === 'number' && (r.userInteractionCount | 0) === cur) {
          try {
            off && off();
          } catch { }
          resolve({ uic: r.userInteractionCount | 0, ok: true });
        }
      });
      // Send with counters set
      ws.setCounters({ userInteractionCount: cur });
      const pos = api.state.positions.map((p) => [p.x, p.y, p.z]);
      ws.userInteraction({ positions: pos });
      setTimeout(() => {
        try {
          off && off();
        } catch { }
        resolve({ uic: -1, ok: false });
      }, 8000);
    });
  });

  expect(echoed && echoed.ok).toBeTruthy();
  if (echoed && echoed.ok) expect(echoed.uic).toBe(987654321);
});
