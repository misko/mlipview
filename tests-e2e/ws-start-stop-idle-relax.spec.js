import { test, expect } from '@playwright/test';

test('autoMD → stop → idle drag → relax → stop (WS)', async ({ page, baseURL }) => {
  test.setTimeout(60_000);
  // Ensure WS base and disable autoMD via query flag to control start explicitly
  await page.addInitScript(() => {
    window.__MLIPVIEW_SERVER = 'http://127.0.0.1:8000';
    window.__MLIPVIEW_TEST_MODE = true;
  });
  // Health precheck: ensure we're talking to the right server and CUDA is available
  const health = await page.request.get('http://127.0.0.1:8000/serve/health');
  if (!health.ok()) test.skip(true, 'Backend health unavailable');
  const h = await health.json();
  if (!h || String(h.device||'').toLowerCase() !== 'cuda') test.skip(true, 'Backend not on CUDA');
  await page.goto(`${baseURL || ''}/index.html?mol=molecules/water.xyz`);
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 45000 });

  // Connect WS and ensure an initial idle energy frame so UMA is warmed up
  await page.evaluate(async ()=>{
    try {
      const ws = window.__fairchem_ws__ || window.__WS_API__;
      await ws.ensureConnected();
      try { await window.viewerApi?.baselineEnergy?.(); } catch {}
      try { await ws.waitForEnergy({ timeoutMs: 20000 }); } catch {}
    } catch(e){ /* swallow errors in page to avoid throwing Event */ }
  });

  // Start MD and collect >= 12 frames
  const mdFrames = await page.evaluate(async ()=>{
    const ws = window.__fairchem_ws__;
    const N = (window.viewerApi?.state?.positions||[]).length|0;
    let count = 0;
    const off = ws.onResult((r)=>{ try { if (Array.isArray(r.positions) && r.positions.length===N) count++; } catch{} });
    try {
      window.viewerApi.startMDContinuous({ steps: 100, temperature: 1200 });
      const t0 = Date.now();
      while (Date.now() - t0 < 25000 && count < 12) await new Promise(r=>setTimeout(r,100));
    } finally { off && off(); }
    return count;
  });
  expect(mdFrames).toBeGreaterThanOrEqual(12);

  // Stop MD; expect a stop indicator (message or flag)
  const sawStopMd = await page.evaluate(async ()=>{
    const ws = window.__fairchem_ws__;
    let saw = false; const off = ws.onResult((r)=>{ try { if (r && (r.message === 'SIMULATION_STOPPED' || r.simulationStopped === true)) saw = true; } catch{} });
    try {
      window.viewerApi.stopSimulation();
      const t0 = Date.now();
      while (Date.now() - t0 < 5000 && !saw) await new Promise(r=>setTimeout(r,50));
    } finally { off && off(); }
    return saw;
  });
  expect(sawStopMd).toBeTruthy();

  // Drag one atom to generate idle USER_INTERACTION updates; expect >=12 idle results (no positions)
  const idleResults = await page.evaluate(async ()=>{
    const ws = window.__fairchem_ws__;
    const N = (window.viewerApi?.state?.positions||[]).length|0;
    let count = 0; const off = ws.onResult((r)=>{ try { if (r && typeof r.energy === 'number' && (!r.positions || r.positions.length !== N)) count++; } catch{} });
    try {
      const idx = 1; // hydrogen
      const start = window.viewerApi.state.positions[idx];
      // select the atom and perform a small drag oscillation via manipulation API
      try { window.viewerApi.selection.clickAtom(idx); } catch {}
      // Begin drag on a plane perpendicular to X axis through current point
      const planePoint = { x:start.x, y:start.y, z:start.z };
      const planeNormal = { x:1, y:0, z:0 };
      const intersector0 = () => ({ x: planePoint.x, y: planePoint.y, z: planePoint.z });
      window.viewerApi.manipulation.beginDrag(intersector0, { planePoint, planeNormal, source: 'desktop' });
      // simulate small oscillations along X by updating the intersector result
      for (let i=0;i<20;i++){
        const dx = (i%2 ? 1 : -1) * 0.01;
        const p = { x: start.x + dx, y: start.y, z: start.z };
        const intersector = () => ({ x: p.x, y: p.y, z: p.z });
        window.viewerApi.manipulation.updateDrag(intersector);
        // deterministically trigger an idle compute frame and ack it
        try { await window.viewerApi.requestSimpleCalculateNow(); } catch {}
        await new Promise(r=>setTimeout(r, 50));
      }
      window.viewerApi.manipulation.endDrag();
      const t0 = Date.now();
      while (Date.now() - t0 < 5000 && count < 12) await new Promise(r=>setTimeout(r,50));
    } finally { off && off(); }
    return count;
  });
  expect(idleResults).toBeGreaterThanOrEqual(12);

  // Start RELAX and collect >=12 frames
  const relaxFrames = await page.evaluate(async ()=>{
    const ws = window.__fairchem_ws__;
    let count = 0; const off = ws.onResult((r)=>{ try { if (Array.isArray(r.positions)) count++; } catch{} });
    try {
      await window.viewerApi.startRelaxContinuous({ maxSteps: 200 });
      const t0 = Date.now();
      while (Date.now() - t0 < 15000 && count < 12) await new Promise(r=>setTimeout(r,100));
    } finally { off && off(); }
    return count;
  });
  expect(relaxFrames).toBeGreaterThanOrEqual(12);

  // Stop RELAX; expect stop indicator
  const sawStopRelax = await page.evaluate(async ()=>{
    const ws = window.__fairchem_ws__;
    let saw = false; const off = ws.onResult((r)=>{ try { if (r && (r.message === 'SIMULATION_STOPPED' || r.simulationStopped === true)) saw = true; } catch{} });
    try {
      window.viewerApi.stopSimulation();
      const t0 = Date.now();
      while (Date.now() - t0 < 6000 && !saw) await new Promise(r=>setTimeout(r,50));
    } finally { off && off(); }
    return saw;
  });
  expect(sawStopRelax).toBeTruthy();
});
