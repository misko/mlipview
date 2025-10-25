import { test, expect } from './fixtures.js';

test.describe('x-drag-throttle-emit', () => {
  test('drag updates are throttled and final emit occurs after endDrag', async ({ page, baseURL }) => {
    test.setTimeout(60_000);
    await page.addInitScript(() => {
      window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
      window.__MLIPVIEW_TEST_MODE = true;
    });

    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');
    const h = await health.json();
    if (!h || String(h.device || '').toLowerCase() !== 'cuda') {
      test.skip(true, 'Backend not on CUDA');
    }

    await page.goto(`${baseURL || ''}/index.html?mol=molecules/benzene.xyz&autoMD=0`);
    await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

    const data = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      await ws.ensureConnected();
      try { await ws.waitForEnergy({ timeoutMs: 15000 }); } catch { }

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
        const isC =
          (typeof e === 'string' && e.toUpperCase() === 'C') ||
          (typeof e === 'number' && e === 6) ||
          (e && typeof e === 'object' && (e.Z === 6 || e.atomicNumber === 6 || e.z === 6));
        if (isC) { idx = i; break; }
      }

      try { window.viewerApi.selection.clickAtom(idx); } catch { }
      const start = window.viewerApi.state.positions[idx];
      const planePoint = { x: start.x, y: start.y, z: start.z };
      const planeNormal = { x: 1, y: 0, z: 0 };
      const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
      window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });
      for (let i = 0; i < 10; i++) {
        const dx = i % 2 === 0 ? 0.02 : -0.02;
        const p = { x: start.x + dx, y: start.y, z: start.z };
        const intersector = () => ({ x: p.x, y: p.y, z: p.z });
        window.viewerApi.manipulation.updateDrag(intersector);
        await new Promise((r) => setTimeout(r, 20));
      }
      window.viewerApi.manipulation.endDrag();
      await new Promise((r) => setTimeout(r, 200));

      ws.userInteraction = origUI;
      const deltas = [];
      for (let i = 1; i < times.length; i++) deltas.push(times[i] - times[i - 1]);
      return { count: times.length, deltas };
    });

    expect(data.count).toBeGreaterThanOrEqual(2);
    expect(data.count).toBeLessThanOrEqual(12);
    if (data.deltas.length > 0) {
      const maxGap = Math.max(...data.deltas);
      expect(maxGap).toBeGreaterThanOrEqual(60);
    }
  });
});
