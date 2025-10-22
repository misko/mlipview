// Playwright e2e: verify simulation frames snapshot UIC from start and simStep increments
// Uses real WS backend started by global-setup.

import { test, expect } from './fixtures.js';

test.describe('WS protocol: simulation counters', () => {
  test('sim frames echo UIC snapshot and increment simStep', async ({ page, baseURL }) => {
    test.setTimeout(45_000);
    await page.goto(`${baseURL || ''}/index.html?autoMD=0`);
    await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
      timeout: 45000,
    });

    const res = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      const api = window.viewerApi;
      await ws.ensureConnected();

      // Set counters before starting simulation
      const UIC = 777;
      const STEP0 = 100;
      ws.setCounters({ userInteractionCount: UIC, simStep: STEP0 });

      // Prime server state: send USER_INTERACTION with positions carrying counters
      try {
        const pos = api.state.positions.map((p) => [p.x, p.y, p.z]);
        ws.userInteraction({ positions: pos });
        // Wait briefly to allow idle compute and counters update
        await new Promise((r) => setTimeout(r, 250));
      } catch { }

      return await new Promise((resolve) => {
        const off = ws.onResult((r) => {
          try {
            if (r && r.seq) ws.ack(r.seq | 0); // keep server window advancing
          } catch { }
          // Look for a simulation frame that matches our UIC snapshot
          if (r && Array.isArray(r.positions) && (r.userInteractionCount | 0) === UIC) {
            const out = { uic: r.userInteractionCount | 0, simStep: r.simStep | 0 };
            try {
              off && off();
            } catch { }
            resolve(out);
          }
        });
        // Start a cheap MD simulation using LJ to avoid UMA dependency
        ws.startSimulation({
          type: 'md',
          params: { calculator: 'uma', temperature: 300, timestep_fs: 1.0, friction: 0.01 },
        });
        setTimeout(() => {
          try {
            off && off();
          } catch { }
          resolve(null);
        }, 35000);
      });
    });

    expect(res).toBeTruthy();
    expect(res.uic).toBe(777);
    // Allow either STEP0 or STEP0+1 depending on when increment occurs relative to send
    expect([100, 101]).toContain(res.simStep);
  });
});
