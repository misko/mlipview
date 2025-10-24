// tests-e2e/ws-backpressure-ack.spec.js
// Purpose: End-to-end validation of WS backpressure and ACK flow (UMA-only).
// Strategy: block ACK frames at the transport layer (WebSocket.send) until we
// observe WAITING_FOR_ACK; then immediately unblock and send a single ACK.
//
// Run example:
//   WS_DEBUG=1 WS_LOG_CALLS=1 UMA_GEOM_DEBUG=1 \
//   npm run -s test:e2e -- tests-e2e/ws-backpressure-ack.spec.js

import { test, expect } from './fixtures.js';

test.describe('WS protocol: backpressure and ACK', () => {
  test('WAITING_FOR_ACK appears and then clears after ack', async ({ page, loadViewerPage }) => {
    test.setTimeout(45_000);

    const ackDropperInit = () => {
      const proto = WebSocket.prototype;
      if (!Object.prototype.hasOwnProperty.call(proto, '_realSend')) {
        Object.defineProperty(proto, '_realSend', {
          value: proto.send,
          configurable: true,
        });
      }

      window.__BLOCK_ACKS__ = true;
      window.__MLIPVIEW_DEBUG_API = true;

      proto.send = function (data) {
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
    };

    await loadViewerPage({
      query: { autoMD: 0, wsDebug: 1 },
      server: 'http://localhost:8000',
      configOverrides: { minStepIntervalMs: 1, mdFriction: 0.02 },
      extraInit: ackDropperInit,
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
    }

    expect(outcome.waitingSeen).toBeTruthy();
    expect(outcome.cleared).toBeTruthy();
    //expect(consoleErrors.join('\n')).toBe('');
  });
});
