import { test, expect } from './fixtures.js';

test('timeline dock visibility and overlay styling', async ({ page, loadViewerPage }) => {
  test.setTimeout(120_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 150, temperature: 1200 });
  });

  await page.waitForFunction((target) => {
    const stats = window.viewerApi.timeline.bufferStats();
    return stats.size >= target;
  }, 30, { timeout: 40_000 });

  await page.evaluate(() => window.viewerApi.stopSimulation());

  const panel = page.locator('[data-testid="timeline-panel"]');
  await panel.waitFor({ state: 'attached' });
  const initialState = await page.evaluate(() => window.viewerApi.timeline.getState());
  expect(initialState.visible).toBe(false);

  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();
  await expect(panel).toBeVisible();

  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.evaluate((node) => {
    node.value = '-8';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(() => window.viewerApi.timeline.getState().active, null, { timeout: 10_000 });

  const overlayFilter = await page.locator('[data-testid="timeline-overlay"]').evaluate((el) => {
    const styles = getComputedStyle(el);
    return styles.backdropFilter || styles.webkitBackdropFilter || el.style.backdropFilter || '';
  });
  expect(overlayFilter === '' || overlayFilter.toLowerCase() === 'none').toBeTruthy();

  await page.click('[data-testid="timeline-live"]');
  await page.waitForFunction(() => !window.viewerApi.timeline.getState().active, null, { timeout: 20_000 });
});
