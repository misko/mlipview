import { test, expect } from './fixtures.js';

test('Drag updates are throttled to ~100ms and final emit occurs on endDrag', async ({ page, baseURL }) => {
  test.setTimeout(60_000);
  await page.addInitScript(() => {
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
    window.__MLIPVIEW_TEST_MODE = true;
  });

  // Backend health precheck
  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');

  await page.goto(`${baseURL || ''}/index.html?mol=molecules/benzene.xyz&autoMD=0`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

  const data = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    await ws.ensureConnected();
    try { await ws.waitForEnergy({ timeoutMs: 15000 }); } catch { }
    // Wrap userInteraction to timestamp calls
    const times = [];
    const origUI = ws.userInteraction.bind(ws);
    ws.userInteraction = (arg) => {
      times.push(performance.now());
      return origUI(arg);
    };

    const els = window.viewerApi?.state?.elements || [];
    let idx = 0;
    for (let i = 0; i < els.length; i++) {
      const e = els[i];
      const isC = (typeof e === 'string' && e.toUpperCase() === 'C') ||
        (typeof e === 'number' && e === 6) ||
        (e && typeof e === 'object' && (e.Z === 6 || e.atomicNumber === 6 || e.z === 6));
      if (isC) { idx = i; break; }
    }
    // Begin drag and spam updates faster than 100ms to exercise throttle
    try { window.viewerApi.selection.clickAtom(idx); } catch { }
    const start = window.viewerApi.state.positions[idx];
    const planePoint = { x: start.x, y: start.y, z: start.z };
    const planeNormal = { x: 1, y: 0, z: 0 };
    const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
    window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });
    for (let i = 0; i < 10; i++) {
      const dx = (i % 2 === 0 ? 0.02 : -0.02);
      const p = { x: start.x + dx, y: start.y, z: start.z };
      const intersector = () => ({ x: p.x, y: p.y, z: p.z });
      window.viewerApi.manipulation.updateDrag(intersector);
      await new Promise(r => setTimeout(r, 20)); // 20ms interval, below throttle window
    }
    window.viewerApi.manipulation.endDrag();
    // Allow a moment for final emit
    await new Promise(r => setTimeout(r, 200));
    ws.userInteraction = origUI;
    // Compute inter-arrival stats
    const deltas = [];
    for (let i = 1; i < times.length; i++) deltas.push(times[i] - times[i - 1]);
    return { count: times.length, deltas };
  });

  // Expect we did not emit for each of the 10 updates; throttle should reduce count roughly <= ~ceil(duration/100ms)+1
  expect(data.count).toBeGreaterThanOrEqual(2);
  expect(data.count).toBeLessThanOrEqual(7);
  // Most deltas should be >= ~80ms (allow some jitter)
  const largeGaps = data.deltas.filter((d) => d >= 80);
  expect(largeGaps.length).toBeGreaterThanOrEqual(Math.max(1, data.deltas.length - 2));
});
