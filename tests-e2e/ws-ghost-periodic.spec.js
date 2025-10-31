import { test, expect } from './fixtures.js';

async function ensureGhostsEnabled(page) {
  await page.evaluate(() => {
    const api = window.viewerApi;
    const { state } = api;
    if (typeof state.toggleCellVisibilityEnhanced === 'function' && !state.showCell) {
      state.toggleCellVisibilityEnhanced();
    }
    if (typeof state.toggleGhostCells === 'function' && !state.showGhostCells) {
      state.toggleGhostCells();
    }
    state.markCellChanged?.();
    api.view.rebuildGhosts?.();
  });
}

async function ghostSummary(page) {
  return await page.evaluate(() => {
    const groups = window.viewerApi.view?._internals?.ghostBondGroups;
    if (!groups) return { size: 0, allSoft: false };
    const entries = Array.from(groups.values());
    const allSoft = entries.every((entry) => entry.master === entry.baseGroup?.masters?.soft);
    return { size: entries.length, allSoft };
  });
}

async function waitForTimelineFrames(page, target) {
  await page.waitForFunction(
    (min) => window.viewerApi.timeline.bufferStats().size >= min,
    target,
    { timeout: 45_000 }
  );
}

async function selectHistoricalOffset(page) {
  return await page.evaluate(() => {
    const offsets = window.viewerApi.timeline.getOffsets();
    if (!offsets || !offsets.length) return null;
    const candidates = offsets.filter((o) => o < -1);
    if (!candidates.length) return offsets[offsets.length - 1];
    return candidates[Math.floor(candidates.length / 2)];
  });
}

test.describe('ghost periodic meshes', () => {
  test('ghost bonds reuse soft masters during live and timeline playback', async ({ page, loadViewerPage }) => {
    test.setTimeout(180_000);

    await loadViewerPage({
      query: { mol: 'molecules/benzene.xyz', autoMD: 0 },
      disableAutoMd: true,
    });

    await ensureGhostsEnabled(page);
    const liveSnapshot = await ghostSummary(page);
    expect(liveSnapshot.size).toBeGreaterThan(0);
    expect(liveSnapshot.allSoft).toBe(true);

    await page.evaluate(() => {
      window.viewerApi.startMDContinuous({ steps: 220, temperature: 1400 });
    });
    await waitForTimelineFrames(page, 30);

    const offset = await selectHistoricalOffset(page);
    expect(offset).not.toBeNull();

    await page.evaluate((value) => window.viewerApi.timeline.select(value), offset);
    await page.waitForFunction(() => window.viewerApi.timeline.getState().active, null, { timeout: 12_000 });

    await page.evaluate(() => window.viewerApi.view.rebuildGhosts?.());

    const timelineSnapshot = await ghostSummary(page);
    expect(timelineSnapshot.size).toBeGreaterThan(0);
    expect(timelineSnapshot.allSoft).toBe(true);

    await page.evaluate(() => window.viewerApi.timeline.live());
    await page.waitForFunction(() => !window.viewerApi.timeline.getState().active, null, { timeout: 12_000 });
    await page.evaluate(() => window.viewerApi.stopSimulation());
  });
});
