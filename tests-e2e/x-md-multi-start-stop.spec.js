// Regression harness for websocket-driven MD runs:
// 1. Repeatedly start/stop continuous MD and collect every incoming frame to prove
//    sequence numbers stay monotonic with no duplicate deliveries.
// 2. Programmatically trigger the MD start action three times in quick succession
//    (bypassing the UI guard) and assert via debugStreamListenerStats that only one
//    ws.onResult listener is ever attached and later detached.
// The viewer exposes the stats hook in test mode so we can verify the wiring without
// poking private internals.
import { test, expect } from './fixtures.js';

const BACKEND_URL = 'http://127.0.0.1:8000';

test('MD start/stop cycles do not duplicate decoded frames', async ({
  page,
  loadViewerPage,
  ensureWsReady,
}) => {
  test.setTimeout(120_000);

  const health = await page.request.get(`${BACKEND_URL}/serve/health`);
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const healthInfo = await health.json();
  if (!healthInfo || String(healthInfo.device || '').toLowerCase() !== 'cuda')
    test.skip(true, 'Backend not on CUDA');

  await loadViewerPage({
    query: { mol: 'molecules/water.xyz', autoMD: 0 },
    server: BACKEND_URL,
    testMode: false,
    disableAutoMd: true,
  });

  expect(await ensureWsReady()).toBeTruthy();

  const result = await page.evaluate(async () => {
    const wait = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
    const ws = window.__fairchem_ws__;
    if (!ws || typeof ws.onResult !== 'function') throw new Error('WS client unavailable');

    const metrics = () => {
      try {
        return window.viewerApi?.getMetrics?.() || {};
      } catch {
        return {};
      }
    };

    const results = [];

    for (let cycle = 0; cycle < 3; cycle++) {
      const record = { frames: 0, duplicates: [], nonMonotonic: [] };
      const seenSeq = new Set();
      let lastSeq = -Infinity;

      const off = ws.onResult((r) => {
        try {
          if (!r || typeof r.seq !== 'number') return;
          if (!Array.isArray(r.positions) || !r.positions.length) return;
          const seq = Number(r.seq);
          record.frames++;
          if (seenSeq.has(seq)) record.duplicates.push(seq);
          else seenSeq.add(seq);
          if (seq <= lastSeq) record.nonMonotonic.push({ seq, lastSeq });
          lastSeq = seq;
        } catch {
          /* ignore */
        }
      });

      try {
        const startRes = await window.viewerApi.startMDContinuous({
          steps: 160,
          temperature: 1200,
        });
        if (startRes?.disabled) throw new Error('MD loop disabled');

        const startDeadline = Date.now() + 8000;
        while (Date.now() < startDeadline) {
          const running = metrics().running;
          if (running) break;
          await wait(100);
        }

        const collectDeadline = Date.now() + 8000;
        while (Date.now() < collectDeadline && record.frames < 12) {
          await wait(100);
        }

        try {
          await window.viewerApi.stopSimulation();
        } catch {
          /* ignore stop errors */
        }

        const stopDeadline = Date.now() + 8000;
        while (Date.now() < stopDeadline) {
          if (!metrics().running) break;
          await wait(100);
        }
      } finally {
        off && off();
      }

      results.push(record);
      await wait(300);
    }

    let multiStart = null;
    try {
      const debugStatsFn = window.viewerApi?.debugStreamListenerStats;
      if (typeof debugStatsFn === 'function') {
        debugStatsFn({ reset: true });
        await window.viewerApi.stopSimulation();
        await wait(300);

        const startResults = [];
        for (let i = 0; i < 3; i++) {
          startResults.push(await window.viewerApi.startMDContinuous({
            steps: 160,
            temperature: 900,
          }));
          await wait(30);
        }

        const startDeadline = Date.now() + 6000;
        while (Date.now() < startDeadline) {
          if (metrics().running) break;
          await wait(100);
        }

        await wait(500);
        await window.viewerApi.stopSimulation();

        const stats = debugStatsFn();
        multiStart = { startResults, stats };
      }
    } catch {
      /* optional multi-start diagnostics */
    }

    return { cycles: results, multiStart };
  });

  expect(result.cycles.length).toBe(3);
  for (const [index, cycle] of result.cycles.entries()) {
    expect(cycle.frames, `cycle ${index} produced frames`).toBeGreaterThanOrEqual(8);
    expect(cycle.duplicates, `cycle ${index} duplicate seqs`).toEqual([]);
    expect(cycle.nonMonotonic, `cycle ${index} non-monotonic seqs`).toEqual([]);
  }

  expect(result.multiStart, 'multi-start diagnostics').not.toBeNull();
  const { startResults, stats } = result.multiStart;
  expect(startResults.length).toBe(3);
  expect(startResults[0]?.streaming).toBeTruthy();
  expect(startResults[1]?.ignored).toBeTruthy();
  expect(startResults[2]?.ignored).toBeTruthy();
  expect(stats?.md?.attach).toBe(1);
  expect(stats?.md?.detach).toBe(1);
  expect(stats?.md?.maxActive).toBe(1);
  expect(stats?.md?.active).toBe(0);
});
