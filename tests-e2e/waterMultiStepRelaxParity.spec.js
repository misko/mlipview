import { test, expect } from './fixtures.js';
import fs from 'fs';

test.setTimeout(60000);

// Configurable tolerance (eV) per-step and overall RMS
const MAX_ABS_PER_STEP = 0.01; // single step deviations should be tiny
const MAX_RMS = 0.005; // cumulative RMS tolerance
const STEPS = 20; // number of relax steps to perform

function loadReference() {
  const text = fs.readFileSync('fairchem_local_server/water_fairchem_bfgs_trace.json', 'utf8');
  const json = JSON.parse(text);
  // Use first STEPS+1 energies (initial + STEPS steps) if available
  return json.energies.slice(0, STEPS + 1);
}

async function ensureWater(page) {
  await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true);
  await page.waitForSelector('#moleculeSelect');
  await page.selectOption('#moleculeSelect', 'molecules/water.xyz');
  // Wait until water coordinates loaded (3 atoms)
  await page.waitForFunction(() => (window.viewerApi?.state?.positions?.length || 0) === 3);
  // Force baseline energy (sync remote fetch) and retry until finite; rely on same-origin proxy (baseURL 4000 -> FastAPI via relative path)
  await page.evaluate(async () => {
    window.__FORCE_SYNC_REMOTE__ = true;
    window.__MLIPVIEW_SERVER = 'http://localhost:8000';
    try {
      await window.viewerApi?.baselineEnergy?.();
    } catch (e) {
      console.log('baselineEnergy error', e);
    }
  });
  // Fallback polling: manually invoke ff.computeForces sync if energy not set yet
  await page.waitForFunction(
    () => {
      const have =
        typeof window.viewerApi?.state?.dynamics?.energy === 'number' &&
        isFinite(window.viewerApi.state.dynamics.energy);
      if (!have) {
        try {
          window.viewerApi?.ff?.computeForces?.({ sync: true });
        } catch {}
      }
      return (
        typeof window.viewerApi?.state?.dynamics?.energy === 'number' &&
        isFinite(window.viewerApi.state.dynamics?.energy)
      );
    },
    null,
    { timeout: 20000 }
  );
  // Confirm geometry atomic numbers
  const geom = await page.evaluate(() => window.__dumpCurrentAtoms?.());
  expect(geom.atomic_numbers).toEqual([1, 1, 8]);
}

test('multi-step relax parity (20 steps) vs reference BFGS trace', async ({ page }) => {
  const ref = loadReference();
  expect(ref.length).toBeGreaterThanOrEqual(STEPS + 1);
  await page.addInitScript(() => {
    window.__FORCE_SYNC_REMOTE__ = true;
  });
  await page.goto('/');
  await ensureWater(page);
  // Collect energies: initial + each step after clicking relax button
  const energies = [];
  const initialE = await page.evaluate(() => window.viewerApi.state.dynamics.energy);
  energies.push(initialE);
  // Perform STEPS relax steps
  for (let i = 0; i < STEPS; i++) {
    await page.click('#btnRelax');
    // wait a tiny bit for energy update to propagate
    await page.waitForTimeout(120);
    const E = await page.evaluate(() => window.viewerApi.state.dynamics.energy);
    energies.push(E);
  }
  // Basic sanity: energies should not diverge wildly
  for (const e of energies) {
    expect(e).toBeLessThan(-1000);
    expect(e).toBeGreaterThan(-3000);
  }
  // Compare per-step deltas to reference
  // ref[0] corresponds to initial reference energy, energies[0] current initial
  const absDiffs = energies.map((e, i) => Math.abs(e - ref[i]));
  absDiffs.forEach((d, i) => {
    // Allow slightly looser tolerance for step 0 because model initialization can shift a few 1e-3
    const tol = i === 0 ? MAX_ABS_PER_STEP * 2 : MAX_ABS_PER_STEP;
    expect(d).toBeLessThanOrEqual(tol);
  });
  const rms = Math.sqrt(absDiffs.reduce((a, b) => a + b * b, 0) / absDiffs.length);
  expect(rms).toBeLessThanOrEqual(MAX_RMS);
  // Log summary for debugging
  console.log(
    'RELAX20_PARITY',
    JSON.stringify({
      initialDiff: absDiffs[0],
      maxDiff: Math.max(...absDiffs),
      rms,
      energies,
      ref: ref.slice(0, energies.length),
    })
  );
  // Persist artifact
  await page.context()._fsWrite?.(); // no-op placeholder if custom hook exists
  try {
    const artifact = {
      steps: STEPS,
      maxAbs: Math.max(...absDiffs),
      rms,
      perStepDiffs: absDiffs,
      energies,
      reference: ref.slice(0, energies.length),
      timestamp: new Date().toISOString(),
    };
    // Use page.evaluate to expose data then save via injected Node fs (Playwright test context is Node, so we can write directly here).
    const fsMod = await import('fs');
    if (!fsMod.existsSync('test-results')) fsMod.mkdirSync('test-results');
    fsMod.writeFileSync(
      'test-results/waterMultiStepRelaxParity.json',
      JSON.stringify(artifact, null, 2)
    );
  } catch (e) {
    console.warn('Failed to write parity artifact', e);
  }
});
