// Purpose: End-to-end verification that idle computes (no simulation running)
// produce frames with energy (and possibly forces) but omit positions as per
// protobuf_migration.md protocol guidance.
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
          } catch { }
        });
        setTimeout(() => {
          try {
            off && off();
          } catch { }
          resolve(out);
        }, timeout);
      });
    },
    { timeout }
  );
}

test('idle frames omit positions and carry energy', async ({ page, baseURL }) => {
  test.setTimeout(45000);
  // Stream full browser console to test output and collect errors for assertions
  const consoleErrors = [];
  page.on('console', (msg) => {
    const text = msg.text();
    if (msg.type() === 'error') consoleErrors.push(text);
    // Prefix with [browser:<type>] so it's easy to scan in CI logs
    // eslint-disable-next-line no-console
    console.log(`[browser:${msg.type()}] ${text}`);
  });
  page.on('pageerror', (err) => {
    const text = (err && (err.message || String(err))) || 'unknown pageerror';
    consoleErrors.push(text);
    // eslint-disable-next-line no-console
    console.log(`[pageerror] ${text}`);
  });
  page.on('requestfailed', (req) => {
    const failure = (req.failure && req.failure()) || {};
    // eslint-disable-next-line no-console
    console.log(
      `[requestfailed] ${req.method()} ${req.url()} -> ${failure.errorText || 'failed'}`
    );
  });
  await page.addInitScript(() => {
    window.__MLIPVIEW_TEST_MODE = false;
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  });
  await page.goto(`${baseURL}/index.html?autoMD=0&wsDebug=1&debug=1`);

  await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
    timeout: 45000,
  });
  // Ensure WS is connected before triggering the idle compute
  await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    if (ws && typeof ws.ensureConnected === 'function') {
      await ws.ensureConnected();
    }
  });
  // Trigger an idle compute, but FIRST attach a listener to avoid race
  const idle = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    const api = window.viewerApi;
    if (!ws || !api) return null;
    await ws.ensureConnected();
    return await new Promise((resolve) => {
      const off = ws.onResult((r) => {
        if (r && typeof r.energy === 'number') {
          try { off && off(); } catch { }
          resolve(r);
        }
      });
      // Now request the idle compute
      api.requestSimpleCalculateNow();
      setTimeout(() => {
        try { off && off(); } catch { }
        resolve(null);
      }, 10000);
    });
  });
  expect(idle).toBeTruthy();
  // Idle frames should not contain positions per protocol
  if (idle) {
    const hasPositions = Array.isArray(idle.positions) && idle.positions.length > 0;
    expect(hasPositions).toBeFalsy();
  }
  // Optional: ensure no protobuf/WS init errors in console
  const errBlob = consoleErrors.join('\n');
  expect(errBlob).not.toMatch(/protobuf-missing|Failed to parse message/);
});
