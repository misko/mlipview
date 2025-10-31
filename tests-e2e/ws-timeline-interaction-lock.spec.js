import { test, expect } from './fixtures.js';

test('timeline mode locks interactions until live resumes', async ({ page, loadViewerPage }) => {
  test.setTimeout(120_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  // Warm up MD stream to populate buffer
  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 200, temperature: 1500 });
  });

  await page.waitForFunction((target) => {
    const stats = window.viewerApi?.timeline?.bufferStats?.();
    return stats && stats.size >= target;
  }, 40, { timeout: 45_000 });

  await page.evaluate(() => window.viewerApi.stopSimulation());

  // Select an atom so manipulation APIs have context
  await page.evaluate(() => window.viewerApi.selection.clickAtom(0));

  // Enter timeline via slider interaction
  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();
  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.evaluate((node) => {
    node.value = '-10';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    return !!status && !!status.active;
  }, null, { timeout: 10_000 });

  const beginDuringTimeline = await page.evaluate(() => window.viewerApi.manipulation.beginDrag(
    () => ({ x: 0, y: 0, z: 0 }),
    { planePoint: { x: 0, y: 0, z: 0 }, planeNormal: { x: 1, y: 0, z: 0 } },
  ));
  expect(beginDuringTimeline).toBe(false);

  await page.click('[data-testid="timeline-live"]');
  await page.waitForFunction(() => {
    const timeline = window.viewerApi.timeline;
    const status =
      (timeline && typeof timeline.getStatus === 'function' && timeline.getStatus()) ||
      (timeline && typeof timeline.getState === 'function' && timeline.getState()) ||
      null;
    return !!status && !status.active;
  }, null, { timeout: 20_000 });

  const beginAfterLive = await page.evaluate(() => window.viewerApi.manipulation.beginDrag(
    () => ({ x: 0, y: 0, z: 0 }),
    { planePoint: { x: 0, y: 0, z: 0 }, planeNormal: { x: 1, y: 0, z: 0 } },
  ));
  expect(beginAfterLive).toBe(true);

  await page.evaluate(() => window.viewerApi.manipulation.endDrag?.());
});
