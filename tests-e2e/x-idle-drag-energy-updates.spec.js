import { test, expect } from './fixtures.js';

test.describe('x-idle-drag-energy-updates', () => {
  test('idle drag emits USER_INTERACTION with positions and refreshes forces/energy', async ({ page, baseURL }) => {
    test.setTimeout(60_000);
    await page.addInitScript(() => {
      window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
    });

    const health = await page.request.get('http://127.0.0.1:8000/serve/health');
    if (!health.ok()) test.skip(true, 'Backend health unavailable');
    const h = await health.json();
    if (!h || String(h.device || '').toLowerCase() !== 'cuda') {
      test.skip(true, 'Backend not on CUDA');
    }

    await page.goto(`${baseURL || ''}/index.html?mol=molecules/benzene.xyz&autoMD=0`);
    await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

    await page.evaluate(async () => {
      const ws = window.__fairchem_ws__ || window.__WS_API__;
      await ws.ensureConnected();
      try { await ws.waitForEnergy({ timeoutMs: 15000 }); } catch { }
      try { await window.viewerApi?.baselineEnergy?.(); } catch { }
    });

    const results = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      const N = (window.viewerApi?.state?.positions || []).length | 0;
      const els = window.viewerApi?.state?.elements || [];
      let idx = 0;
      for (let i = 0; i < els.length; i++) {
        const e = els[i];
        const isC =
          (typeof e === 'string' && e.toUpperCase() === 'C') ||
          (typeof e === 'number' && e === 6) ||
          (e && typeof e === 'object' && (e.Z === 6 || e.atomicNumber === 6 || e.z === 6));
        const isH =
          (typeof e === 'string' && e.toUpperCase() === 'H') ||
          (typeof e === 'number' && e === 1) ||
          (e && typeof e === 'object' && (e.Z === 1 || e.atomicNumber === 1 || e.z === 1));
        if (isC) { idx = i; break; }
        if (!isH) idx = i;
      }

      let sentUI = 0;
      let sentUIWithPos = 0;
      let energyFrames = 0;
      let forcesUpdated = false;
      const energyLen0 = window.viewerApi?.debugEnergySeriesLength?.() || 0;

      const origUI = ws.userInteraction.bind(ws);
      ws.userInteraction = (arg) => {
        try {
          sentUI++;
          if (arg && Array.isArray(arg.positions)) sentUIWithPos++;
        } catch {}
        return origUI(arg);
      };
      const off = ws.onResult((r) => {
        try {
          if (r && typeof r.energy === 'number') energyFrames++;
        } catch {}
      });

      try { window.viewerApi.selection.clickAtom(idx); } catch {}
      const start = window.viewerApi.state.positions[idx];
      const planePoint = { x: start.x, y: start.y, z: start.z };
      const planeNormal = { x: 1, y: 0, z: 0 };
      const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
      window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });
      for (let i = 0; i < 20; i++) {
        const dx = i % 4 === 0 ? 0.02 : i % 4 === 1 ? -0.02 : i % 4 === 2 ? 0.015 : -0.015;
        const p = { x: start.x + dx, y: start.y, z: start.z };
        const intersector = () => ({ x: p.x, y: p.y, z: p.z });
        window.viewerApi.manipulation.updateDrag(intersector);
        try { await window.viewerApi.requestSimpleCalculateNow(); } catch {}
        await new Promise((r) => setTimeout(r, 30));
        try {
          const f = window.viewerApi?.state?.forces;
          if (Array.isArray(f) && f.length === N) forcesUpdated = true;
        } catch {}
      }
      window.viewerApi.manipulation.endDrag();

      const t0 = Date.now();
      while (Date.now() - t0 < 1500) {
        try {
          const f = window.viewerApi?.state?.forces;
          if (Array.isArray(f) && f.length === N) { forcesUpdated = true; break; }
        } catch {}
        await new Promise((r) => setTimeout(r, 50));
      }

      const energyLen1 = window.viewerApi?.debugEnergySeriesLength?.() || energyLen0;
      off && off();
      ws.userInteraction = origUI;
      return {
        sentUI,
        sentUIWithPos,
        energyFrames,
        forcesUpdated,
        energyTicksDelta: energyLen1 - energyLen0,
      };
    });

    expect(results.sentUI).toBeGreaterThanOrEqual(6);
    expect(results.sentUIWithPos).toBeGreaterThanOrEqual(6);
    expect(results.energyFrames).toBeGreaterThanOrEqual(3);
    expect(results.forcesUpdated).toBeTruthy();
    expect(results.energyTicksDelta).toBeGreaterThanOrEqual(2);
  });
});

