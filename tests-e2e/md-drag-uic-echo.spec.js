import { test, expect } from '@playwright/test';

test('MD: after drag, server frames echo non-zero UIC', async ({ page, baseURL }) => {
  test.setTimeout(60_000);
  await page.addInitScript(() => {
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  });

  // Backend health precheck
  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');

  await page.goto(`${baseURL || ''}/index.html?mol=molecules/benzene.xyz&autoMD=0`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

  // Start MD
  await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    await ws.ensureConnected();
    try { await ws.waitForEnergy({ timeoutMs: 15000 }); } catch { }
    await window.viewerApi.startMDContinuous({ steps: 1000, temperature: 1000 });
  });

  const results = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    const els = window.viewerApi?.state?.elements || [];
    let idx = 0;
    for (let i = 0; i < els.length; i++) {
      const e = els[i];
      const isC = (typeof e === 'string' && e.toUpperCase() === 'C') ||
        (typeof e === 'number' && e === 6) ||
        (e && typeof e === 'object' && (e.Z === 6 || e.atomicNumber === 6 || e.z === 6));
      if (isC) { idx = i; break; }
    }

    let uicBefore = [];
    let uicAfter = [];
    let off = ws.onResult((r) => {
      try {
        if (typeof r.userInteractionCount === 'number') {
          // Collect samples categorized by time window
          if (!window.__DRAG_DONE__) uicBefore.push(r.userInteractionCount);
          else uicAfter.push(r.userInteractionCount);
        }
      } catch { }
    });

    // Give MD a moment to run
    await new Promise(r => setTimeout(r, 500));

    // Perform a small drag while MD runs; this should send USER_INTERACTION with current counters
    try { window.viewerApi.selection.clickAtom(idx); } catch { }
    const start = window.viewerApi.state.positions[idx];
    const planePoint = { x: start.x, y: start.y, z: start.z };
    const planeNormal = { x: 1, y: 0, z: 0 };
    const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
    window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });
    for (let i = 0; i < 6; i++) {
      const dx = (i % 2 === 0 ? 0.02 : -0.02);
      const p = { x: start.x + dx, y: start.y, z: start.z };
      const intersector = () => ({ x: p.x, y: p.y, z: p.z });
      window.viewerApi.manipulation.updateDrag(intersector);
      await new Promise(r => setTimeout(r, 25));
    }
    window.viewerApi.manipulation.endDrag();
    window.__DRAG_DONE__ = true;

    // Wait for some frames after drag
    await new Promise(r => setTimeout(r, 800));
    off && off();
    // Compute stats
    const beforeMax = uicBefore.length ? Math.max(...uicBefore) : 0;
    const afterMax = uicAfter.length ? Math.max(...uicAfter) : 0;
    return { beforeMax, afterMax, beforeSamples: uicBefore.slice(0, 5), afterSamples: uicAfter.slice(0, 5) };
  });

  // Expect that after the drag, server frames begin echoing a non-zero UIC (or a strictly larger value)
  expect(results.afterMax).toBeGreaterThanOrEqual(1);
});
