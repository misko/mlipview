import { test, expect } from './fixtures.js';

test.describe('x-relax-drag-latch', () => {
  test('drag during relax stays latched on next frame', async ({ page, baseURL }) => {
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
      try { await ws.waitForEnergy({ timeoutMs: 15000 }); } catch {}
    });

    const kept = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      const N = (window.viewerApi?.state?.positions || []).length | 0;
      const els = window.viewerApi?.state?.elements || [];
      const carbonIdx = Array.isArray(els)
        ? els.findIndex((e) => {
            if (typeof e === 'string') return e.toUpperCase() === 'C';
            if (typeof e === 'number') return e === 6;
            if (e && typeof e === 'object') return e.Z === 6 || e.atomicNumber === 6 || e.z === 6;
            return false;
          })
        : -1;
      const idx = carbonIdx >= 0 ? carbonIdx : 0;
      let frameCount = 0;
      let draggedTo = null;
      const samples = [];
      let resolved = false;
      let resolveOuter;
      const done = new Promise((r) => (resolveOuter = r));
      const off = ws.onResult((r) => {
        try {
          if (!r || !Array.isArray(r.positions) || r.positions.length !== N) return;
          frameCount++;
          if (frameCount === 12) {
            const start = window.viewerApi.state.positions[idx];
            try { window.viewerApi.selection.clickAtom(idx); } catch {}
            const planePoint = { x: start.x, y: start.y, z: start.z };
            const planeNormal = { x: 1, y: 0, z: 0 };
            const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
            window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });
            const p = { x: start.x + 0.02, y: start.y, z: start.z };
            const intersector = () => ({ x: p.x, y: p.y, z: p.z });
            window.viewerApi.manipulation.updateDrag(intersector);
            window.viewerApi.manipulation.endDrag();
            draggedTo = p;
          } else if (frameCount > 12 && draggedTo) {
            const cur = window.viewerApi.state.positions[idx];
            samples.push({ x: cur.x, y: cur.y, z: cur.z });
            if (!resolved && samples.length >= 3) {
              resolved = samples.some((curPos) => {
                const dx = Math.abs(curPos.x - draggedTo.x);
                const dy = Math.abs(curPos.y - draggedTo.y);
                const dz = Math.abs(curPos.z - draggedTo.z);
                return dx < 1e-2 && dy < 1e-2 && dz < 1e-2;
              });
              if (resolved) resolveOuter(true);
            }
          }
        } catch {}
      });
      try {
        await window.viewerApi.startRelaxContinuous({ maxSteps: 500 });
      } catch {}
      setTimeout(() => {
        if (!resolved) resolveOuter(false);
      }, 25000);
      const ok = await done;
      off && off();
      return ok;
    });

    expect(kept).toBeTruthy();
  });
});
