import { test, expect, request } from '@playwright/test';
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

async function buildReferenceTrace(){
  // Use backend /relax sequential one-step calls to mimic UI steps and store energies
  const { symbols, coords } = loadROYXYZ();
  const atomic_numbers = symbols.map(symbolToZ);
  const refEnergies = [];
  // First: single-point baseline energy via /simple_calculate
  const singleResp = await fetch('http://localhost:8000/simple_calculate', {
    method:'POST', headers:{'Content-Type':'application/json'},
    body: JSON.stringify({ atomic_numbers, coordinates: coords, properties:['energy','forces'], calculator:'uma' })
  });
  if(!singleResp.ok) throw new Error('Reference baseline failed '+singleResp.status);
  const singleJson = await singleResp.json();
  refEnergies.push(singleJson.results.energy);
  let currentCoords = coords.map(r=> [...r]);
  for(let i=0;i<STEPS;i++){
    const resp = await fetch('http://localhost:8000/relax', { method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify({ atomic_numbers, coordinates: currentCoords, steps:1, calculator:'uma' }) });
    if(!resp.ok) throw new Error('Reference relax step failed '+resp.status);
    const js = await resp.json();
    refEnergies.push(js.final_energy);
    currentCoords = js.positions; // update for next step
  }
  return refEnergies;
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

test('ROY multi-step relax parity (20 steps) vs dynamic backend reference', async ({ page }) => {
  // Build reference energies (initial + steps)
  const ref = await buildReferenceTrace();
  expect(ref.length).toBe(STEPS+1);
  await page.goto('/');
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
