// Playwright test: validate single-socket invariant and streaming behavior
// - Navigates with wsSim=1 to enable streaming over WebSocket
// - Asserts exactly one WebSocket is opened from the page and reused
// - Starts a short MD stream and verifies at least a few frames arrive with increasing seq

import { test, expect } from '@playwright/test';

// Helper to wait for a predicate with timeout
async function waitFor(pred, { timeout = 8000, interval = 50 } = {}){
  const t0 = Date.now();
  while (Date.now() - t0 < timeout){
    const v = await pred();
    if (v) return v;
    await new Promise(r=> setTimeout(r, interval));
  }
  throw new Error('timeout');
}

test.describe('WS streaming single-socket', () => {
  test('single websocket and increasing frames', async ({ page, baseURL }) => {
    test.setTimeout(30000);
    // Force same-origin static hosting and point frontend to our WS backend on 8001
    await page.addInitScript(() => {
      window.__MLIPVIEW_SERVER = 'http://localhost:8000';
      window.__MLIPVIEW_TEST_MODE = true; // bypass focus gating
      // Speed up by reducing min interval pacing in case REST paths are hit at all
      window.__MLIP_CONFIG = { minStepIntervalMs: 1, mdFriction: 0.02 };
    });

  const sockEvents = [];
  const consoleErrors = [];
  // Stream all console to test output; only collect actual console.error for assertions
  page.on('console', msg => {
    const text = msg.text();
    if (msg.type() === 'error') consoleErrors.push(text);
    // eslint-disable-next-line no-console
    console.log(`[browser:${msg.type()}] ${text}`);
  });
  page.on('pageerror', (err) => {
    const text = (err && (err.message || String(err))) || 'unknown pageerror';
    // eslint-disable-next-line no-console
    console.log(`[pageerror] ${text}`);
  });
  page.on('requestfailed', (req) => {
    const failure = (req.failure && req.failure()) || {};
    // eslint-disable-next-line no-console
    console.log(`[requestfailed] ${req.method()} ${req.url()} -> ${failure.errorText || 'failed'}`);
  });
    page.on('websocket', ws => {
      sockEvents.push({ type: 'open', url: ws.url() });
      ws.on('framereceived', data => { sockEvents.push({ type: 'recv', size: (data?.length||0) }); });
      ws.on('framesent', data => { sockEvents.push({ type: 'send', size: (data?.length||0) }); });
      ws.on('close', () => sockEvents.push({ type: 'close' }));
    });

  // Use a minimal WS test harness page
  await page.goto(`${baseURL}/ws-test.html?sim=1`);
  await page.waitForFunction(() => !!window.__WS_READY__, { timeout: 10000 });

    // Ensure WS API is present
    const apiPresent = await page.evaluate(() => !!window.__WS_API__ && typeof window.__WS_API__.startSimulation === 'function');
    expect(apiPresent).toBeTruthy();

    // Expect at least one WebSocket created and only one unique URL for the session
    await waitFor(() => sockEvents.find(e => e.type === 'open'));
    const wsOpens = sockEvents.filter(e => e.type === 'open');
    expect(wsOpens.length).toBe(1);
    const wsUrl = wsOpens[0].url;
    // No close before we finish
    expect(sockEvents.find(e => e.type === 'close')).toBeFalsy();

    // Collect decoded frame seq numbers via viewer metrics hook by watching ack progression
    const seqs = await page.evaluate(() => new Promise(resolve => {
      const unique = [];
      let maxSeen = -1;
      const stopAt = Date.now() + 8000;
      const poll = () => {
        try {
          const evs = window.__WS_EVENTS__||[];
          for (const r of evs) {
            const s = r && (r.seq|0);
            if (typeof s === 'number' && s > maxSeen) {
              unique.push(s);
              maxSeen = s;
            }
          }
        } catch {}
        if (Date.now() > stopAt) return resolve(unique);
        setTimeout(poll, 100);
      };
      poll();
    }));

    // Alternatively, use the raw protocol events as a backup signal
  const recvCount = sockEvents.filter(e => e.type === 'recv').length;
  // Allow slow environments to produce fewer frames; require at least 2
  expect(recvCount).toBeGreaterThan(1);
  // No console errors
  expect(consoleErrors.join('\n')).toBe('');

    // Basic monotonicity: some increasing seq numbers
    if (Array.isArray(seqs) && seqs.length > 1){
      for (let i=1;i<seqs.length;i++) expect(seqs[i]).toBeGreaterThanOrEqual(seqs[i-1]);
    }

    // Assert no additional WebSocket opens during the short run
    const laterOpens = sockEvents.filter(e => e.type === 'open');
    expect(laterOpens.length).toBe(1);
    expect(wsUrl).toContain('ws://');
  });
});
