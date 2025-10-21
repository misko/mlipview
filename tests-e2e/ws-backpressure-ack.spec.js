// tests-e2e/ws-backpressure-ack.spec.js
// Purpose: End-to-end validation of WS backpressure and ACK flow (UMA-only).
// Strategy: block ACK frames at the transport layer (WebSocket.send) until we
// observe WAITING_FOR_ACK; then immediately unblock and send a single ACK.
//
// Run example:
//   WS_DEBUG=1 WS_LOG_CALLS=1 UMA_GEOM_DEBUG=1 \
//   npm run -s test:e2e -- tests-e2e/ws-backpressure-ack.spec.js

import { test, expect } from '@playwright/test';

test.describe('WS protocol: backpressure and ACK', () => {
  test('WAITING_FOR_ACK appears and then clears after ack', async ({ page, baseURL }) => {
    test.setTimeout(45_000);

    // --- Transport-level ACK dropper BEFORE the app loads ---
    await page.addInitScript(() => {
      const _realSend = WebSocket.prototype.send;
      Object.defineProperty(WebSocket.prototype, '_realSend', {
        value: _realSend,
        configurable: true,
      });

      // When true, heuristic-drops tiny binary frames (typical ACK len 9–11).
      // We’ll flip this to false once we decide to allow ACKs.
      window.__BLOCK_ACKS__ = true;

      WebSocket.prototype.send = function (data) {
        if (window.__BLOCK_ACKS__) {
          try {
            const isBinary = typeof data !== 'string';
            const len =
              (data && typeof data.byteLength === 'number' && data.byteLength) ||
              (data && typeof data.size === 'number' && data.size) ||
              0;

            if (isBinary && len > 0 && len <= 12) {
              console.log(`[test] dropping tiny binary frame len=${len} (likely ACK)`);
              return;
            }
          } catch {
            // if heuristic fails, just send
          }
        }
        return WebSocket.prototype._realSend.call(this, data);
      };
    });

    // --- App flags BEFORE navigation ---
    await page.addInitScript(() => {
      window.__MLIPVIEW_SERVER = 'http://localhost:8000';
      window.__MLIPVIEW_NO_AUTO_MD = true;   // don't auto-start MD
      window.__MLIPVIEW_TEST_MODE = true;    // bypass focus gating
      window.__MLIP_CONFIG = { minStepIntervalMs: 1, mdFriction: 0.02 };
      window.__MLIPVIEW_DEBUG_API = true;    // extra logs if client supports
    });

    // --- Pipe browser logs to your terminal ---
    const consoleErrors = [];
    page.on('console', (msg) => {
      const text = msg.text();
      if (msg.type() === 'error') consoleErrors.push(text);
      console.log(`[browser:${msg.type()}] ${text}`);
    });
    page.on('pageerror', (err) => {
      const text = (err && (err.message || String(err))) || 'unknown pageerror';
      console.log(`[pageerror] ${text}`);
      consoleErrors.push(text);
    });
    page.on('requestfailed', (req) => {
      const failure = (req.failure && req.failure()) || {};
      console.log(
        `[requestfailed] ${req.method()} ${req.url()} -> ${failure.errorText || 'failed'}`
      );
    });

    // --- Load app and wait ready ---
    await page.goto(`${baseURL || ''}/index.html?autoMD=0&debug=1&wsDebug=1`);
    await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
      timeout: 45_000,
    });

    // --- Core test ---
    const outcome = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      await ws.ensureConnected();

      // Belt & suspenders: disable any high-level auto-acks if present.
      if (ws.sendAck) ws.sendAck = () => { };
      if (ws._sendAck) ws._sendAck = () => { };
      try { if ('autoAck' in ws) ws.autoAck = false; } catch { }

      let waitingSeen = false;
      let waitingSeq = 0;
      let cleared = false;
      let firstSimSeq = 0;
      let lastSeq = 0;
      let ackSent = false;

      return await new Promise((resolve) => {
        const off = ws.onResult((r) => {
          try {
            const seq = (r && r.seq) | 0;
            lastSeq = seq;

            // First sim frame usually carries positions array
            if (!firstSimSeq && Array.isArray(r?.positions)) {
              firstSimSeq = seq;
              if (window.__MLIPVIEW_DEBUG_API) console.log('[test] firstSimSeq', firstSimSeq);
            }

            if (r && r.message === 'WAITING_FOR_ACK') {
              // Backpressure detected
              if (!waitingSeen) {
                waitingSeen = true;
                waitingSeq = seq;
                if (window.__MLIPVIEW_DEBUG_API) console.log('[test] WAITING_FOR_ACK at', waitingSeq);

                // Immediately unblock ACKs and send one to clear the stall.
                // Note: server stops advancing seq during stall, so waiting to pass
                // "firstSimSeq + N" can deadlock. Ack now.
                window.__BLOCK_ACKS__ = false;
                if (typeof ws.ack === 'function') {
                  ws.ack(seq);
                } else if (typeof ws.ping === 'function') {
                  ws.ping({ ack: seq });
                }
                ackSent = true;
                if (window.__MLIPVIEW_DEBUG_API) console.log('[test] sent ACK', seq);
              }
            } else if (waitingSeen && ackSent && seq > waitingSeq) {
              // Got a non-waiting message after sending ACK → cleared
              cleared = true;
              try { off && off(); } catch { }
              resolve({ waitingSeen, cleared, waitingSeq, lastSeq, firstSimSeq });
            }
          } catch { }
        });

        // Start UMA MD streaming
        ws.startSimulation({
          type: 'md',
          params: { calculator: 'uma', temperature: 300, timestep_fs: 1.0, friction: 0.01 },
        });

        // Safety timeout
        setTimeout(() => {
          try { off && off(); } catch { }
          resolve({ waitingSeen, cleared, waitingSeq, lastSeq, firstSimSeq, timeout: true });
        }, 35_000);
      });
    });

    if (!outcome.waitingSeen || !outcome.cleared) {
      console.log('[e2e] outcome', outcome);
      console.log('[e2e] consoleErrors:\n' + consoleErrors.join('\n'));
    }

    expect(outcome.waitingSeen).toBeTruthy();
    expect(outcome.cleared).toBeTruthy();
    //expect(consoleErrors.join('\n')).toBe('');
  });
});
