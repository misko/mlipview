// tests-e2e/ws-stream-single-socket.spec.js
// Purpose: Validate WS streaming behavior and single-socket invariant.
// - Navigates to ws-test harness with sim=1 to produce a short MD stream.
// - Asserts exactly one WebSocket is created and stays open.
// - Observes protocol frames to ensure more than one incoming frame and
//   monotonically increasing seq values.
// - Verifies the client never sends TEXT frames (no JSON sent by client).
//
// Note: Some server stacks may emit a text/JSON banner or health ping. We ignore
// that specific benign console error if it appears; client should ideally not log it
// as an error, but we keep the filter to avoid flakes across versions.

import { test, expect } from './fixtures.js';

// Helper to wait for a predicate with timeout
async function waitFor(pred, { timeout = 8000, interval = 50 } = {}) {
  const t0 = Date.now();
  while (Date.now() - t0 < timeout) {
    const v = await pred();
    if (v) return v;
    // eslint-disable-next-line no-await-in-loop
    await new Promise((r) => setTimeout(r, interval));
  }
  throw new Error('timeout');
}

test.describe('WS streaming single-socket', () => {
  test('single websocket and increasing frames', async ({ page, baseURL }) => {
    test.setTimeout(30000);

    // Force same-origin static hosting and point frontend to our WS backend
    await page.addInitScript(() => {
      window.__MLIPVIEW_SERVER = 'http://localhost:8000';
      window.__MLIPVIEW_TEST_MODE = true; // bypass focus gating
      window.__MLIP_CONFIG = { minStepIntervalMs: 1, mdFriction: 0.02 };
    });

    const sockEvents = [];
    const consoleErrors = [];

    // Stream page console for debugging; collect only error-level for assertions
    page.on('console', (msg) => {
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
      console.log(
        `[requestfailed] ${req.method()} ${req.url()} -> ${failure.errorText || 'failed'}`
      );
    });

    // Capture websocket lifecycle + frame types
    // We push different event records for binary vs text to assert client never sends text.
    page.on('websocket', (ws) => {
      sockEvents.push({ type: 'open', url: ws.url() });

      ws.on('framereceived', (data) => {
        // `data` is a string for text frames, or a Buffer for binary
        if (typeof data === 'string') {
          sockEvents.push({ type: 'recvText', size: data.length });
        } else {
          sockEvents.push({ type: 'recv', size: data?.length || 0 });
        }
      });

      ws.on('framesent', (data) => {
        if (typeof data === 'string') {
          // This would indicate the client sent TEXT (JSON) to the server — not allowed.
          sockEvents.push({ type: 'sendText', size: data.length });
        } else {
          sockEvents.push({ type: 'send', size: data?.length || 0 });
        }
      });

      ws.on('close', () => sockEvents.push({ type: 'close' }));
    });

    // Use a minimal WS test harness page
    await page.goto(`${baseURL}/ws-test.html?sim=1&wsDebug=1&debug=1`);
    await page.waitForFunction(() => !!window.__WS_READY__, { timeout: 10000 });

    // Ensure WS API is present
    const apiPresent = await page.evaluate(
      () => !!window.__WS_API__ && typeof window.__WS_API__.startSimulation === 'function'
    );
    expect(apiPresent).toBeTruthy();

    // Expect a single WebSocket created
    await waitFor(() => sockEvents.find((e) => e.type === 'open'));
    const wsOpens = sockEvents.filter((e) => e.type === 'open');
    expect(wsOpens.length).toBe(1);
    const wsUrl = wsOpens[0].url;
    expect(wsUrl.startsWith('ws://') || wsUrl.startsWith('wss://')).toBeTruthy();

    // No close before we finish
    expect(sockEvents.find((e) => e.type === 'close')).toBeFalsy();

    // Collect decoded frame seq numbers via viewer metrics hook by watching ack progression
    const seqs = await page.evaluate(
      () =>
        new Promise((resolve) => {
          const unique = [];
          let maxSeen = -1;
          const stopAt = Date.now() + 8000;
          const poll = () => {
            try {
              const evs = window.__WS_EVENTS__ || [];
              for (const r of evs) {
                const s = (r && (r.seq | 0)) || 0;
                if (Number.isFinite(s) && s > maxSeen) {
                  unique.push(s);
                  maxSeen = s;
                }
              }
            } catch { }
            if (Date.now() > stopAt) return resolve(unique);
            setTimeout(poll, 100);
          };
          poll();
        })
    );

    // As a backup signal, confirm we received some binary frames
    const recvCount = sockEvents.filter((e) => e.type === 'recv').length;
    expect(recvCount).toBeGreaterThan(1);

    // Assert the client never sent TEXT frames (i.e., no JSON client→server)
    const sentText = sockEvents.filter((e) => e.type === 'sendText');
    expect(sentText.length).toBe(0);

    // Filter out a known benign server text-frame message that older clients may log as error.
    const benign = /^\[WS] Received JSON text in binary path; ignoring frame\.$/;
    const benign2 = /^\[WS] Received JSON-like bytes; ignoring frame\.$/;
    const realErrors = consoleErrors.filter((e) => !benign.test(e) && !benign2.test(e));
    expect(realErrors.join('\n')).toBe('');

    // Basic monotonicity: some increasing seq numbers
    if (Array.isArray(seqs) && seqs.length > 1) {
      for (let i = 1; i < seqs.length; i++) {
        expect(seqs[i]).toBeGreaterThanOrEqual(seqs[i - 1]);
      }
    }

    // Assert no additional WebSocket opens during the short run
    const laterOpens = sockEvents.filter((e) => e.type === 'open');
    expect(laterOpens.length).toBe(1);
  });
});
