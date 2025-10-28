// Timeline rewind/resume acceptance test
import { test, expect } from './fixtures.js';

test('timeline rewind preserves frame identity and resumes live stream', async ({ page, loadViewerPage }) => {
  test.setTimeout(120_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  // Kick off an MD run long enough to populate the frame buffer.
  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 200, temperature: 1200 });
  });

  await page.waitForFunction((target) => {
    const stats = window.viewerApi?.timeline?.bufferStats?.();
    return stats && stats.size >= target;
  }, 60, { timeout: 45_000 });

  const initialSize = await page.evaluate(() => window.viewerApi.timeline.bufferStats().size);

  // Reveal the timeline dock.
  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();
  await expect(page.locator('[data-testid="timeline-panel"]')).toBeVisible();

  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.waitFor({ state: 'visible' });

  // Select frame -15.
  await slider.evaluate((node) => {
    node.value = '-15';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(() => {
    const state = window.viewerApi.timeline.getState();
    return state.active && state.offset === -15;
  }, null, { timeout: 10_000 });

  const sigFirst = await page.evaluate(() => window.viewerApi.timeline.getSignature(-15));
  expect(sigFirst).toBeTruthy();

  // Move to frame -16.
  await slider.evaluate((node) => {
    node.value = '-16';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(() => window.viewerApi.timeline.getState().offset === -16, null, { timeout: 10_000 });

  // Return to frame -15 and confirm identical payload.
  await slider.evaluate((node) => {
    node.value = '-15';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(() => window.viewerApi.timeline.getState().offset === -15, null, { timeout: 10_000 });
  const sigSecond = await page.evaluate(() => window.viewerApi.timeline.getSignature(-15));
  expect(sigSecond).toBe(sigFirst);

  // Play forward from the selected frame and ensure we land on the latest frame.
  await page.click('[data-testid="timeline-play"]');
  await page.waitForFunction(() => {
    const state = window.viewerApi.timeline.getState();
    return !state.active && state.mode === 'live' && state.offset === -1;
  }, null, { timeout: 30_000 });

  await page.waitForFunction(() => window.viewerApi.getMetrics().running === 'md', null, { timeout: 20_000 });

  const targetSize = initialSize + 30;
  await page.waitForFunction((expected) => {
    const stats = window.viewerApi.timeline.bufferStats();
    return stats.size >= expected;
  }, targetSize, { timeout: 60_000 });

  await page.evaluate(() => window.viewerApi.stopSimulation());
});
