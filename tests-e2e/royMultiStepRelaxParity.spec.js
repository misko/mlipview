import { test, expect } from '@playwright/test';
import fs from 'fs';

test.setTimeout(60000);

const STEPS = 20;
const MAX_ABS_PER_STEP = 0.01;
const MAX_RMS = 0.005;

// Utility to fetch ROY xyz from repo
function loadROYXYZ(){
  const txt = fs.readFileSync('public/molecules/roy.xyz','utf8');
  const lines = txt.trim().split(/\r?\n/).filter(l=>l.trim().length>0);
  const n = parseInt(lines[0].trim(),10);
  const atoms = lines.slice(2,2+n).map(l=> l.trim().split(/\s+/));
  const symbols = atoms.map(a=> a[0]);
  const coords = atoms.map(a=> [parseFloat(a[1]),parseFloat(a[2]),parseFloat(a[3])]);
  return { symbols, coords };
}

// Minimal symbol->Z map (reuse subset from frontend)
const SYMBOL_TO_Z = {H:1, He:2, C:6, N:7, O:8, S:16};
function symbolToZ(s){ return SYMBOL_TO_Z[s] || 0; }

async function buildReferenceTrace(page){
  // Build reference energies using the WS client directly (no REST)
  const { symbols, coords } = loadROYXYZ();
  const atomic_numbers = symbols.map(symbolToZ);
  // Navigate to the minimal WS harness to ensure WS client is available
  // and not influenced by the full viewer's state.
  await page.addInitScript(() => { window.__MLIPVIEW_SERVER = 'http://localhost:8000'; window.__MLIPVIEW_TEST_MODE = true; });
  await page.goto('/ws-test.html?sim=0');
  await page.waitForFunction(() => !!window.__WS_READY__, null, { timeout: 10000 });
  // Perform INIT_SYSTEM and WS simple calculate, then run a streaming relax and collect N frames.
  const energies = await page.evaluate(async ({ atomic_numbers, coords, steps }) => {
    const ws = window.__WS_API__ || window.__fairchem_ws__;
    if (!ws) throw new Error('WS API not present on harness page');
    await ws.ensureConnected();
    ws.initSystem({ atomic_numbers, positions: coords });
    // Baseline energy via WS simple calculate
    const single = await ws.requestSimpleCalculate();
    const out = [single.energy];
    // Start a streaming relax and collect energies from subsequent frames
    window.__WS_EVENTS__ = [];
    const off = ws.onResult((r)=>{ try { if (r && r.seq) ws.ack(r.seq|0); } catch{} });
    ws.startSimulation({ type: 'relax', params: { calculator: 'uma', max_step: 0.01 } });
    const t0 = Date.now();
    while (out.length < (steps + 1)) {
      const evs = window.__WS_EVENTS__ || [];
      for (const r of evs) {
        if (r && typeof r.energy === 'number') {
          // append only new energies; avoid duplicates by checking length vs seen
          if (out.length === 0 || r.energy !== out[out.length - 1]) out.push(r.energy);
        }
      }
      window.__WS_EVENTS__ = [];
      if (Date.now() - t0 > 30000) break; // safety cap 30s
      await new Promise(r=> setTimeout(r, 50));
    }
    try { off && off(); } catch {}
    try { ws.stopSimulation(); } catch {}
    return out.slice(0, steps + 1);
  }, { atomic_numbers, coords, steps: STEPS });
  return energies;
}

async function ensureROY(page){
  await page.waitForFunction(()=>window.__MLIP_DEFAULT_LOADED===true);
  await page.waitForSelector('#moleculeSelect');
  await page.selectOption('#moleculeSelect','molecules/roy.xyz');
  await page.waitForFunction(()=> (window.viewerApi?.state?.positions?.length||0) === 27);
  await page.evaluate(async ()=>{ window.__FORCE_SYNC_REMOTE__=true; window.__MLIPVIEW_SERVER='http://localhost:8000'; try { await window.viewerApi?.baselineEnergy?.(); } catch{} });
  // Poll energy to become finite
  await page.waitForFunction(()=> typeof window.viewerApi?.state?.dynamics?.energy === 'number' && isFinite(window.viewerApi.state.dynamics.energy), null, { timeout:20000 });
}

test('ROY multi-step relax parity (20 steps) vs dynamic backend reference', async ({ page, baseURL }) => {
  // Build reference energies (initial + steps) using WS-only harness
  const ref = await buildReferenceTrace(page);
  expect(ref.length).toBe(STEPS+1);
  // Now load the full viewer UI and perform the same relax steps via UI
  await page.goto(baseURL || '/');
  await ensureROY(page);
  const energies = [];
  energies.push(await page.evaluate(()=> window.viewerApi.state.dynamics.energy));
  for(let i=0;i<STEPS;i++){
    await page.click('#btnRelax');
    await page.waitForTimeout(150);
    energies.push(await page.evaluate(()=> window.viewerApi.state.dynamics.energy));
  }
  // Assertions
  const absDiffs = energies.map((e,i)=> Math.abs(e - ref[i]));
  absDiffs.forEach((d,i)=>{ const tol = i===0 ? MAX_ABS_PER_STEP*2 : MAX_ABS_PER_STEP; expect(d).toBeLessThanOrEqual(tol); });
  const rms = Math.sqrt(absDiffs.reduce((a,b)=>a+b*b,0)/absDiffs.length);
  expect(rms).toBeLessThanOrEqual(MAX_RMS);
  const artifact = { molecule:'ROY', steps:STEPS, maxAbs: Math.max(...absDiffs), rms, perStepDiffs: absDiffs, energies, reference: ref, timestamp: new Date().toISOString() };
  if(!fs.existsSync('test-results')) fs.mkdirSync('test-results');
  fs.writeFileSync('test-results/royMultiStepRelaxParity.json', JSON.stringify(artifact,null,2));
  console.log('ROY_RELAX20_PARITY', JSON.stringify(artifact));
});
