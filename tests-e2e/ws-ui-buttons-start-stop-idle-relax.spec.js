// Purpose: End-to-end UI control test using buttons to start/stop MD and RELAX under WS-only backend.
import { test, expect } from './fixtures.js';

// Uses UI toggles (visible buttons) instead of direct API calls:
// - Auto MD runs on page load (do not enable test mode)
// - Stop via MD toggle and verify stop message/flag
// - Drag (via manipulation API) to generate >=12 idle frames
// - Start Relax via Relax toggle, wait >=12 frames
// - Stop via Relax toggle and verify stop message/flag
test('UI buttons: autoMD → stop → idle drag → relax → stop', async ({
  page,
  loadViewerPage,
}) => {
  test.setTimeout(60_000);

  // Health precheck: backend must be up and on CUDA for UMA
  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');

  await loadViewerPage({
    query: { mol: 'molecules/water.xyz' },
    server: 'http://127.0.0.1:8000',
    testMode: false,
    disableAutoMd: false,
  });

  // Connect WS and warm up one idle energy to ensure UMA is responsive
  await page.evaluate(async () => {
    const ws = window.__fairchem_ws__ || window.__WS_API__;
    await ws.ensureConnected();
    try {
      await ws.waitForEnergy({ timeoutMs: 20000 });
    } catch { }
  });

  // Expect auto MD to run and produce >= 12 frames
  const mdFrames = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    const N = (window.viewerApi?.state?.positions || []).length | 0;
    let count = 0;
    const off = ws.onResult((r) => {
      try {
        if (Array.isArray(r.positions) && r.positions.length === N) count++;
      } catch { }
    });
    try {
      const t0 = Date.now();
      while (Date.now() - t0 < 25000 && count < 12) await new Promise((r) => setTimeout(r, 100));
    } finally {
      off && off();
    }
    return count;
  });
  expect(mdFrames).toBeGreaterThanOrEqual(12);

  // Stop MD via UI toggle and verify stop indicator (use DOM click to avoid visibility issues)
  await page.waitForSelector('#toggleMD', { state: 'attached' });
  const sawStopMd = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    let saw = false;
    const off = ws.onResult((r) => {
      try {
        if (r && (r.message === 'SIMULATION_STOPPED' || r.simulationStopped === true)) saw = true;
      } catch { }
    });
    try {
      const btn = document.getElementById('toggleMD');
      btn?.click();
      const t0 = Date.now();
      while (Date.now() - t0 < 6000 && !saw) await new Promise((r) => setTimeout(r, 50));
    } finally {
      off && off();
    }
    return saw;
  });
  expect(sawStopMd).toBeTruthy();

  // Perform a small drag via manipulation API, requesting deterministic idle compute each step; expect >=12 idle frames
  const idleResults = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    const N = (window.viewerApi?.state?.positions || []).length | 0;
    let count = 0;
    const off = ws.onResult((r) => {
      try {
        if (r && typeof r.energy === 'number' && (!r.positions || r.positions.length !== N))
          count++;
      } catch { }
    });
    try {
      const idx = 1; // hydrogen
      const start = window.viewerApi.state.positions[idx];
      try {
        window.viewerApi.selection.clickAtom(idx);
      } catch { }
      const planePoint = { x: start.x, y: start.y, z: start.z };
      const planeNormal = { x: 1, y: 0, z: 0 };
      const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
      window.viewerApi.manipulation.beginDrag(intersector0, {
        planePoint,
        planeNormal,
        source: 'desktop',
      });
      for (let i = 0; i < 20; i++) {
        const dx = (i % 2 ? 1 : -1) * 0.01;
        const p = { x: start.x + dx, y: start.y, z: start.z };
        const intersector = () => ({ x: p.x, y: p.y, z: p.z });
        window.viewerApi.manipulation.updateDrag(intersector);
        try {
          await window.viewerApi.requestSimpleCalculateNow();
        } catch { }
        await new Promise((r) => setTimeout(r, 50));
      }
      window.viewerApi.manipulation.endDrag();
      const t0 = Date.now();
      while (Date.now() - t0 < 5000 && count < 12) await new Promise((r) => setTimeout(r, 50));
    } finally {
      off && off();
    }
    return count;
  });
  expect(idleResults).toBeGreaterThanOrEqual(12);

  // Start Relax via UI toggle and collect >= 12 frames
  await page.waitForSelector('#toggleRelax', { state: 'attached' });
  await page.evaluate(() => {
    const el = document.getElementById('toggleRelax');
    if (el) el.click();
  });
  const relaxFrames = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    let count = 0;
    const off = ws.onResult((r) => {
      try {
        if (Array.isArray(r.positions)) count++;
      } catch { }
    });
    try {
      const t0 = Date.now();
      while (Date.now() - t0 < 20000 && count < 12) await new Promise((r) => setTimeout(r, 100));
    } finally {
      off && off();
    }
    return count;
  });
  expect(relaxFrames).toBeGreaterThanOrEqual(12);

  // Stop Relax via UI toggle and verify stop indicator
  const sawStopRelax = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    let saw = false;
    const off = ws.onResult((r) => {
      try {
        if (r && (r.message === 'SIMULATION_STOPPED' || r.simulationStopped === true)) saw = true;
      } catch { }
    });
    try {
      const btn = document.getElementById('toggleRelax');
      btn?.click();
      const t0 = Date.now();
      while (Date.now() - t0 < 6000 && !saw) await new Promise((r) => setTimeout(r, 50));
    } finally {
      off && off();
    }
    return saw;
  });
  expect(sawStopRelax).toBeTruthy();
});

test('Temperature indicator only reports during MD runs', async ({
  page,
  loadViewerPage,
  ensureWsReady,
}) => {
  test.setTimeout(60_000);

  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device || '').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');

  await loadViewerPage({
    query: { mol: 'molecules/water.xyz', autoMD: 0 },
    server: 'http://127.0.0.1:8000',
  });

  expect(await ensureWsReady()).toBeTruthy();

  const instTemp = page.locator('#instTemp');
  const placeholderRegex = /T:\s*[—-]\s*K/;

  await expect(instTemp).toHaveText(placeholderRegex);

  await page.evaluate(async () => {
    await window.viewerApi.startMDContinuous({ steps: 500, temperature: 1200 });
  });

  await expect(instTemp).toHaveText(/T:\s*\d+(?:\.\d+)?\s*K/, { timeout: 20_000 });

  await page.evaluate(() => {
    window.viewerApi.stopSimulation();
  });

  await expect(instTemp).toHaveText(placeholderRegex, { timeout: 20_000 });

  const sawRelaxFrames = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    if (!ws) return false;
    await window.viewerApi.startRelaxContinuous({ maxSteps: 120 });
    return await new Promise((resolve) => {
      const off = ws.onResult((r) => {
        if (Array.isArray(r.positions)) {
          try { off && off(); } catch { }
          resolve(true);
        }
      });
      setTimeout(() => {
        try { off && off(); } catch { }
        resolve(false);
      }, 10_000);
    });
  });

  expect(sawRelaxFrames).toBeTruthy();

  await expect(instTemp).toHaveText(placeholderRegex, { timeout: 10_000 });

  await page.evaluate(() => {
    window.viewerApi.stopSimulation();
  });
});
