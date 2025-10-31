import { test, expect } from './fixtures.js';
import {
  waitForTimelineBuffer,
  hoverTimelinePanel,
  selectTimelineOffset,
  waitForTimelineState,
  computeTimelineOffset,
} from './utils/timeline.js';

test('timeline dock visibility and overlay styling', async ({ page, loadViewerPage }) => {
  test.setTimeout(120_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 150, temperature: 1200 });
  });

  await waitForTimelineBuffer(page, 30, { label: 'visibility-buffer' });

  await page.evaluate(() => window.viewerApi.stopSimulation());

  const panel = page.locator('[data-testid="timeline-panel"]');
  await panel.waitFor({ state: 'attached' });
  const initialState = await page.evaluate(() => {
    const timeline = window.viewerApi.timeline;
    return (
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null
    );
  });
  expect(initialState?.visible ?? false).toBe(false);

  await hoverTimelinePanel(page, 'visibility-hover');
  await expect(panel).toBeVisible();

  const offset = await computeTimelineOffset(page, { offset: -8 });
  expect(offset).not.toBeNull();
  await selectTimelineOffset(page, offset, { label: 'visibility-select' });
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
    { label: 'visibility-active' }
  );

  const overlayFilter = await page.locator('[data-testid="timeline-overlay"]').evaluate((el) => {
    const styles = getComputedStyle(el);
    return styles.backdropFilter || styles.webkitBackdropFilter || el.style.backdropFilter || '';
  });
  expect(overlayFilter === '' || overlayFilter.toLowerCase() === 'none').toBeTruthy();

  await page.click('[data-testid="timeline-live"]');
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
    { label: 'visibility-return-live' }
  );
});
