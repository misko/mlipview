import { test, expect } from './fixtures.js';

// Repro for: idle responses have lower UIC than current drag UIC; energy plot should still tick for each response.
test('Idle drag: accept lower-UIC energy frames and tick plot', async ({ page, baseURL }) => {
  test.setTimeout(60_000);
  await page.addInitScript(() => {
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
    // Enable API debug to mirror user report logs if needed
    const q = new URLSearchParams(window.location.search || '');
    if (q.get('debug') == null) {
      const u = new URL(window.location.href);
      u.searchParams.set('debug', '1');
      history.replaceState({}, '', u.toString());
      window.__MLIPVIEW_DEBUG_API = true;
    }
  });

  // Backend health precheck
  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');

  await page.goto(`${baseURL || ''}/index.html?mol=molecules/benzene.xyz&autoMD=0&debug=1`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

  // Warm WS and seed baseline energy
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
    // pick a non-hydrogen if available; prefer carbon
    let idx = 0;
    for (let i = 0; i < els.length; i++) {
      const e = els[i];
      const isC = (typeof e === 'string' && e.toUpperCase() === 'C') ||
        (typeof e === 'number' && e === 6) ||
        (e && typeof e === 'object' && (e.Z === 6 || e.atomicNumber === 6 || e.z === 6));
      const isH = (typeof e === 'string' && e.toUpperCase() === 'H') ||
        (typeof e === 'number' && e === 1) ||
        (e && typeof e === 'object' && (e.Z === 1 || e.atomicNumber === 1 || e.z === 1));
      if (isC) { idx = i; break; }
      if (!isH) idx = i; // fallback
    }

    const energyLen0 = Array.isArray(window.__RELAX_TRACE) ? window.__RELAX_TRACE.length : 0;
    let energyFrames = 0;
    let minServerUIC = Infinity;
    let maxClientUIC = -Infinity;
    const off = ws.onResult((r) => {
      try {
        if (r && typeof r.energy === 'number') energyFrames++;
        if (typeof r.userInteractionCount === 'number') {
          minServerUIC = Math.min(minServerUIC, r.userInteractionCount);
        }
      } catch { }
    });

    // Drag without calling requestSimpleCalculateNow so responses may lag UIC
    try { window.viewerApi.selection.clickAtom(idx); } catch { }
    const start = window.viewerApi.state.positions[idx];
    const planePoint = { x: start.x, y: start.y, z: start.z };
    const planeNormal = { x: 1, y: 0, z: 0 };
    const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
    window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });

    for (let i = 0; i < 16; i++) {
      const dx = (i % 4 === 0 ? 0.02 : i % 4 === 1 ? -0.02 : i % 4 === 2 ? 0.015 : -0.015);
      const p = { x: start.x + dx, y: start.y, z: start.z };
      const intersector = () => ({ x: p.x, y: p.y, z: p.z });
      window.viewerApi.manipulation.updateDrag(intersector);
      await new Promise(r => setTimeout(r, 35));
      try {
        // Track current client UIC for comparison
        const ver = window.viewerApi.getVersionInfo();
        if (ver && typeof ver.userInteractionVersion === 'number') {
          maxClientUIC = Math.max(maxClientUIC, ver.userInteractionVersion);
        }
      } catch { }
    }
    window.viewerApi.manipulation.endDrag();

    // Allow some time for last responses
    await new Promise(r => setTimeout(r, 1000));
    off && off();
    const energyLen1 = Array.isArray(window.__RELAX_TRACE) ? window.__RELAX_TRACE.length : energyLen0;
    return { energyFrames, energyTicksDelta: energyLen1 - energyLen0, minServerUIC, maxClientUIC };
  });

  // We expect several energy frames and at least a couple of plot ticks.
  expect(results.energyFrames).toBeGreaterThanOrEqual(2);
  expect(results.energyTicksDelta).toBeGreaterThanOrEqual(2);
  // In this scenario, it's possible minServerUIC < maxClientUIC; ensure we still plotted.
  expect(results.minServerUIC <= results.maxClientUIC).toBeTruthy();
});
