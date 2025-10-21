// Purpose: End-to-end validation of WS backpressure and ACK flow.
// Ensures that when the client withholds ACKs during a running MD stream,
// the server emits a notice message "WAITING_FOR_ACK" and that once the
// client acknowledges a later seq, the notices stop and streaming resumes.
import { test, expect } from '@playwright/test';

test.describe('WS protocol: backpressure and ACK', () => {
  test('WAITING_FOR_ACK appears and then clears after ack', async ({ page, baseURL }) => {
    test.setTimeout(45_000);
    await page.goto(`${baseURL || ''}/index.html?autoMD=0`);
    await page.waitForFunction(() => !!window.viewerApi && !!window.__MLIP_DEFAULT_LOADED, {
      timeout: 45000,
    });

    const outcome = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      await ws.ensureConnected();
      let waitingSeen = false;
      let seqSeen = 0;
      let firstSimSeq = 0;
      let cleared = false;
      const events = [];
      return await new Promise((resolve) => {
        const off = ws.onResult((r) => {
          try {
            // Intentionally do NOT ack initially to trigger backpressure
            events.push({ seq: r.seq | 0, msg: r.message || null });
            if (r && r.message === 'WAITING_FOR_ACK') {
              waitingSeen = true;
              seqSeen = r.seq | 0;
            }
            if (!firstSimSeq && Array.isArray(r.positions)) firstSimSeq = r.seq | 0;
            if (waitingSeen && !cleared && firstSimSeq && r.seq > firstSimSeq + 15) {
              // Send an ack to clear the backlog
              ws.ack(r.seq | 0);
            }
            if (waitingSeen && !cleared && r.message !== 'WAITING_FOR_ACK' && r.seq > seqSeen) {
              cleared = true;
              try {
                off && off();
              } catch { }
              resolve({ waitingSeen, cleared });
            }
          } catch { }
        });
        ws.startSimulation({
          type: 'md',
          params: { calculator: 'lj', temperature: 300, timestep_fs: 1.0, friction: 0.01 },
        });
        setTimeout(() => {
          try {
            off && off();
          } catch { }
          resolve({ waitingSeen, cleared });
        }, 35000);
      });
    });

    expect(outcome.waitingSeen).toBeTruthy();
    expect(outcome.cleared).toBeTruthy();
  });
});
