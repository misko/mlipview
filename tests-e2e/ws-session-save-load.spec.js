// Verifies that session snapshots capture state, restore geometry, and integrate with reset button.
import { test, expect } from './fixtures.js';

test('session snapshot restores geometry, resumes MD, and reset uses baseline', async ({ page, loadViewerPage }) => {
  test.setTimeout(45_000);
  await loadViewerPage({ query: { mol: 'molecules/benzene.xyz' } });

  const baseline = await page.evaluate(async () => {
    await window.viewerApi.requestSimpleCalculateNow?.();
    const count = window.viewerApi.state.elements.length;
    const firstX = window.viewerApi.state.positions?.[0]?.x || 0;
    const energyLen = window.viewerApi.debugEnergySeriesLength?.() || 0;
    const snapshot = window.viewerApi.session.captureSnapshot({ kind: 'test', label: 'baseline' });
    return { snapshot, count, firstX, energyLen };
  });
  expect(baseline.count).toBeGreaterThan(0);

  const removedCount = await page.evaluate(async () => {
    await window.viewerApi.removeAtomByIndex(0);
    return window.viewerApi.state.elements.length;
  });
  expect(removedCount).toBe(baseline.count - 1);

  await page.evaluate(async (snapshot) => {
    await window.viewerApi.session.loadSnapshot(snapshot);
  }, baseline.snapshot);

  const mdSnapshot = await page.evaluate(async () => {
    const ws = window.__fairchem_ws__;
    let mdFrames = 0;
    const off = ws.onResult((r) => {
      try {
        if (Array.isArray(r.positions)) mdFrames++;
      } catch { }
    });
    try {
      await window.viewerApi.startMDContinuous({ steps: 200, temperature: 600 });
      const t0 = Date.now();
      while (Date.now() - t0 < 12000 && mdFrames < 6) await new Promise((r) => setTimeout(r, 100));
      const snapshot = window.viewerApi.session.captureSnapshot({ kind: 'md', label: 'running-md' });
      return { snapshot, mdFrames };
    } finally {
      off();
    }
  });
  expect(mdSnapshot.mdFrames).toBeGreaterThanOrEqual(2);

  await page.evaluate(() => window.viewerApi.stopSimulation());

  const restored = await page.evaluate(async (snapshot) => {
    await window.viewerApi.session.loadSnapshot(snapshot);
    const timeline = window.viewerApi.timeline;
    const before =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    const ws = window.__fairchem_ws__;
    let frames = 0;
    const off = ws.onResult((r) => {
      try {
        if (Array.isArray(r.positions)) frames++;
      } catch { }
    });
    try {
      // Exit timeline mode; MD should resume automatically after the buffered playback finishes.
      if (timeline && typeof timeline.live === 'function') {
        await timeline.live();
      }
      const t0 = Date.now();
      while (Date.now() - t0 < 12000 && frames < 6) await new Promise((r) => setTimeout(r, 100));
    } finally {
      off();
    }
    return {
      timelineMode: before?.mode ?? 'paused',
      count: window.viewerApi.state.elements.length,
      firstX: window.viewerApi.state.positions?.[0]?.x || 0,
      energyLen: window.viewerApi.debugEnergySeriesLength?.() || 0,
      mdFramesAfterResume: frames,
    };
  }, mdSnapshot.snapshot);

  expect(['paused', 'timeline']).toContain(restored.timelineMode);
  expect(restored.count).toBe(baseline.count);
  expect(restored.energyLen).toBeGreaterThanOrEqual(baseline.energyLen);
  expect(Math.abs(restored.firstX - baseline.firstX)).toBeLessThan(1e-6);
  expect(restored.mdFramesAfterResume).toBeGreaterThanOrEqual(2);

  const mutatedX = await page.evaluate(() => {
    const pos = window.viewerApi.state.positions?.[0];
    if (pos) {
      pos.x += 0.42;
      window.viewerApi.state.markPositionsChanged?.();
      return pos.x;
    }
    return null;
  });
  expect(mutatedX).not.toBeNull();

  await page.evaluate(() => window.viewerApi.session.resetToLastLoad());

  const afterReset = await page.evaluate(() => ({
    count: window.viewerApi.state.elements.length,
    firstX: window.viewerApi.state.positions?.[0]?.x || 0,
  }));
  expect(afterReset.count).toBe(baseline.count);
  expect(Math.abs(afterReset.firstX - baseline.firstX)).toBeLessThan(1e-6);
});
