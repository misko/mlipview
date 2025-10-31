import { test, expect } from './fixtures.js';

test.describe('library query param', () => {
  test('loads SN2 session when ?library=sn2 is present', async ({ page, loadViewerPage }) => {
    test.setTimeout(60_000);

    await loadViewerPage({ query: { library: 'sn2' }, disableAutoMd: false });

    await page.waitForFunction(() => window.__MLIP_DEFAULT_LOADED === true, null, { timeout: 20_000 });

    const info = await page.evaluate(() => {
      const api = window.viewerApi;
      const timelineState = api?.timeline?.getState?.() || {};
      const status = document.getElementById('status')?.textContent || '';
      const stats = api?.timeline?.bufferStats?.() || { size: 0 };
      return {
        atomCount: Array.isArray(api?.state?.positions) ? api.state.positions.length : 0,
        timelineMode: timelineState.mode || null,
        timelineBuffer: stats.size || 0,
        status,
      };
    });

    expect(info.atomCount).toBeGreaterThan(0);
    expect(info.status).toContain('library:sn2');
    expect(info.timelineBuffer).toBeGreaterThan(0);
    expect(typeof info.timelineMode).toBe('string');
  });
});
