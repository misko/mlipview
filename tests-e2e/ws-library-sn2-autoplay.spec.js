import { test, expect } from './fixtures.js';

test.describe('library SN2 playback', () => {
  test('loads SN2 session and resumes MD after playback', async ({ page, loadViewerPage }) => {
    test.setTimeout(60_000);

    await loadViewerPage({ disableAutoMd: false });

    await page.evaluate(async () => {
      const api = window.viewerApi;
      if (!api?.session?.loadFromLibrary) throw new Error('loadFromLibrary unavailable');
      await api.session.loadFromLibrary('sn2');
    });

    await page.waitForFunction(
      () => {
        const state = window.viewerApi?.timeline?.getState?.();
        return state && state.mode === 'playing';
      },
      null,
      { timeout: 15_000 }
    );

    const initialTimeline = await page.evaluate(() => window.viewerApi.timeline.getState());
    expect(initialTimeline.mode).toBe('playing');

    const liveFrames = await page.evaluate(async () => {
      const ws = window.__fairchem_ws__;
      if (!ws?.onResult) return 0;
      let frames = 0;
      const off = ws.onResult((res) => {
        try {
          if (Array.isArray(res.positions)) frames++;
        } catch {}
      });
      try {
        await window.viewerApi.timeline.live();
        const start = Date.now();
        while (Date.now() - start < 12000 && frames < 15) {
          await new Promise((resolve) => setTimeout(resolve, 100));
        }
      } finally {
        off?.();
      }
      return frames;
    });

    expect(liveFrames).toBeGreaterThanOrEqual(15);

    const finalState = await page.evaluate(() => window.__fairchem_ws__?.getState?.() || null);
    expect(finalState?.connected).toBe(true);
  });
});
