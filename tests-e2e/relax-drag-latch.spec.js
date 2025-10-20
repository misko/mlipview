import { test, expect } from '@playwright/test';

test('Relax: drag carbon at 12th frame is latched on next frame', async ({ page, baseURL }) => {
  test.setTimeout(60_000);
  // Ensure we target the right backend
  await page.addInitScript(() => {
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
  });

  // Health precheck: backend must be up and CUDA-enabled
  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');

  await page.goto(`${baseURL || ''}/index.html?mol=molecules/benzene.xyz&autoMD=0`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

  // Warm WS and initial energy
  await page.evaluate(async () => {
    const ws = window.__fairchem_ws__ || window.__WS_API__;
    await ws.ensureConnected();
    try { await ws.waitForEnergy({ timeoutMs: 15000 }); } catch { }
  });

  // Start Relaxation and count frames; after 12 frames, drag a carbon atom and verify the subsequent frame keeps it
  const kept = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    const N = (window.viewerApi?.state?.positions || []).length | 0;
    const els = window.viewerApi?.state?.elements || [];
    const carbonIdx = (Array.isArray(els) ? els.findIndex(e => {
      if (typeof e === 'string') return e.toUpperCase() === 'C';
      if (typeof e === 'number') return e === 6;
      if (e && typeof e === 'object') return (e.Z === 6) || (e.atomicNumber === 6) || (e.z === 6);
      return false;
    }) : -1);
    const idx = carbonIdx >= 0 ? carbonIdx : 0; // fallback to first atom
    let frameCount = 0;
    let afterDragChecked = false;
    let draggedTo = null;
    let resolveOuter;
    const done = new Promise(r => (resolveOuter = r));
    const off = ws.onResult((r) => {
      try {
        if (!r || !Array.isArray(r.positions) || r.positions.length !== N) return;
        frameCount++;
        if (frameCount === 12) {
          // Perform a small desktop drag on the carbon atom
          const start = window.viewerApi.state.positions[idx];
          try { window.viewerApi.selection.clickAtom(idx); } catch { }
          const planePoint = { x: start.x, y: start.y, z: start.z };
          const planeNormal = { x: 1, y: 0, z: 0 };
          const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
          window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'vr' });
          const dx = 0.02; // small offset along x
          const p = { x: start.x + dx, y: start.y, z: start.z };
          const intersector = () => ({ x: p.x, y: p.y, z: p.z });
          window.viewerApi.manipulation.updateDrag(intersector);
          window.viewerApi.manipulation.endDrag();
          draggedTo = p;
        } else if (frameCount > 12 && !afterDragChecked) {
          // On the very next frame after drag, verify the viewer state kept the dragged position
          const cur = window.viewerApi.state.positions[idx];
          const dx = Math.abs(cur.x - draggedTo.x);
          const dy = Math.abs(cur.y - draggedTo.y);
          const dz = Math.abs(cur.z - draggedTo.z);
          const tol = 5e-3; // be lenient to float noise
          afterDragChecked = true;
          resolveOuter(dx < tol && dy < tol && dz < tol);
        }
      } catch { }
    });
    try {
      await window.viewerApi.startRelaxContinuous({ maxSteps: 500 });
    } catch { }
    // In case relax finishes very quickly, also set a timeout safeguard
    setTimeout(() => { if (!afterDragChecked) resolveOuter(false); }, 25000);
    const ok = await done;
    off && off();
    return ok;
  });

  expect(kept).toBeTruthy();
});
