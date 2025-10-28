import { test, expect } from './fixtures.js';

const hoverTimeline = async (page) => {
  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();
  await expect(page.locator('[data-testid="timeline-panel"]')).toBeVisible();
};

const waitForTimelineState = async (page, predicate, opts = {}) => {
  await page.waitForFunction(predicate, opts.arg, { timeout: opts.timeout ?? 20_000 });
};

test('timeline control policy enforcement', async ({ page, loadViewerPage }) => {
  test.setTimeout(180_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 400, temperature: 1500 });
  });

  await waitForTimelineState(page, (target) => {
    const stats = window.viewerApi.timeline.bufferStats();
    return stats.size >= target;
  }, { arg: 60, timeout: 45_000 });

  await hoverTimeline(page);

  const playBtn = page.locator('[data-testid="timeline-play"]');
  const pauseBtn = page.locator('[data-testid="timeline-pause"]');
  const liveBtn = page.locator('[data-testid="timeline-live"]');

  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeEnabled();
  await expect(liveBtn).toBeDisabled();

  await pauseBtn.click();

  await waitForTimelineState(page, () => {
    const state = window.viewerApi.timeline.getState();
    return state.active && state.mode === 'paused' && state.offset === -1;
  });

  await expect(playBtn).toBeEnabled();
  await expect(pauseBtn).toBeDisabled();
  await expect(liveBtn).toBeEnabled();

  await liveBtn.click();

  await waitForTimelineState(page, () => !window.viewerApi.timeline.getState().active, { timeout: 20_000 });

  await hoverTimeline(page);
  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeEnabled();
  await expect(liveBtn).toBeDisabled();

  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.evaluate((node) => {
    node.value = '-12';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await waitForTimelineState(page, () => {
    const state = window.viewerApi.timeline.getState();
    return state.active && state.mode === 'paused' && state.offset === -12;
  });

  await expect(playBtn).toBeEnabled();
  await expect(pauseBtn).toBeDisabled();
  await expect(liveBtn).toBeEnabled();

  await playBtn.click();

  await waitForTimelineState(page, () => window.viewerApi.timeline.getState().mode === 'playing');

  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeEnabled();
  await expect(liveBtn).toBeEnabled();

  await pauseBtn.click();

  await waitForTimelineState(page, () => window.viewerApi.timeline.getState().mode === 'paused');

  await playBtn.click();

  await waitForTimelineState(page, () => window.viewerApi.timeline.getState().mode === 'playing');

  await waitForTimelineState(page, () => {
    const state = window.viewerApi.timeline.getState();
    return state.mode === 'live' && !state.active;
  }, { timeout: 30_000 });

  await hoverTimeline(page);
  await expect(playBtn).toBeDisabled();
  await expect(pauseBtn).toBeEnabled();
  await expect(liveBtn).toBeDisabled();

  await page.evaluate(() => window.viewerApi.stopSimulation());
});
