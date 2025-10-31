import { test, expect } from './fixtures.js';

test('timeline playback surfaces an energy marker', async ({ page, loadViewerPage }) => {
  test.setTimeout(180_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 260, temperature: 1500 });
  });

  await page.waitForFunction(
    (target) => {
      const stats = window.viewerApi.timeline.bufferStats();
      const series = window.viewerApi.debugEnergySeriesLength();
      return stats.size >= target && series >= target;
    },
    40,
    { timeout: 45_000 },
  );

  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();
  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.evaluate((node) => {
    node.value = '-15';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(
    () => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && status.active && status.offset <= -15;
    },
    null,
    { timeout: 10_000 },
  );

  await page.waitForFunction(
    () => !!window.viewerApi.getEnergyMarker(),
    null,
    { timeout: 5_000 },
  );
  const marker = await page.evaluate(() => window.viewerApi.getEnergyMarker());
  expect(marker).toBeTruthy();
  expect(marker.index).toBeGreaterThanOrEqual(0);
  expect(marker.index).toBeLessThan(marker.length);
  expect(typeof marker.energy).toBe('number');

  await page.click('[data-testid="timeline-live"]');
  await page.waitForFunction(
    () => {
      const timeline = window.viewerApi.timeline;
      const status =
        (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
        (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
        null;
      return status && !status.active;
    },
    null,
    { timeout: 20_000 },
  );

  const after = await page.evaluate(() => window.viewerApi.getEnergyMarker());
  expect(after).toBeNull();
});
