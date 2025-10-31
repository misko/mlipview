import { test, expect } from './fixtures.js';
import {
  waitForTimelineBuffer,
  hoverTimelinePanel,
  selectTimelineOffset,
  computeTimelineOffset,
} from './utils/timeline.js';

test.describe('timeline slider selection', () => {
test('single slider change selects requested frame', async ({ page, loadViewerPage, uiMode }) => {
    test.setTimeout(150_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

    await page.evaluate(() => {
      window.viewerApi.startMDContinuous({ steps: 220, temperature: 1600 });
    });

    await waitForTimelineBuffer(page, 40, { timeout: 45_000 });

    await hoverTimelinePanel(page);
    const slider = page.locator('[data-testid="timeline-slider"]');
    await slider.waitFor({ state: 'visible' });

    const targetOffset = await computeTimelineOffset(page, { mode: 'middle' });
    expect(targetOffset).not.toBeNull();

    await page.evaluate(() => {
      window.__changeCount = 0;
      const sliderEl = document.querySelector('[data-testid="timeline-slider"]');
      sliderEl?.addEventListener('change', () => { window.__changeCount++; }, { once: false });
    });

    await selectTimelineOffset(page, targetOffset);

    await page.waitForFunction(
      (value) => {
        const timeline = window.viewerApi.timeline;
        const status =
          (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
          (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
          null;
        return status && status.offset === value && status.mode !== 'live';
      },
      targetOffset,
      { timeout: 10_000 }
    );

    const result = await page.evaluate(() => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return {
        status,
        change: window.__changeCount || 0,
      };
    });
    expect(result.status?.offset).toBe(targetOffset);
    expect(result.status?.mode).toBe('paused');
    const minChanges = uiMode === 'react' ? 0 : 1;
    expect(result.change).toBeGreaterThanOrEqual(minChanges);
  });
});
