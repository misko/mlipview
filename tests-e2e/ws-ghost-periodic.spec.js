import { test, expect } from './fixtures.js';
import {
  waitForTimelineBuffer,
  computeTimelineOffset,
  selectTimelineOffset,
  waitForTimelineState,
} from './utils/timeline.js';

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
    await waitForTimelineBuffer(page, 30, { label: 'ghost-buffer' });

    const offset = await computeTimelineOffset(page, { mode: 'middle' });
    expect(offset).not.toBeNull();

    await selectTimelineOffset(page, offset, { label: 'ghost-select' });
    await waitForTimelineState(
      page,
      () => {
        const status =
          window.viewerApi?.timeline?.getStatus?.() ??
          window.viewerApi?.timeline?.getState?.() ??
          null;
        return !!status && !!status.active;
      },
      null,
      { label: 'ghost-active' }
    );

    await page.evaluate(() => window.viewerApi.view.rebuildGhosts?.());

    const timelineSnapshot = await ghostSummary(page);
    expect(timelineSnapshot.size).toBeGreaterThan(0);
    expect(timelineSnapshot.allSoft).toBe(true);

    await page.evaluate(async () => {
      await (window.viewerApi.timeline.live?.() || Promise.resolve());
    });
    await waitForTimelineState(
      page,
      () => {
        const status =
          window.viewerApi?.timeline?.getStatus?.() ??
          window.viewerApi?.timeline?.getState?.() ??
          null;
        return !!status && !status.active;
      },
      null,
      { label: 'ghost-return-live' }
    );
    await page.evaluate(() => window.viewerApi.stopSimulation());
  });
});
