import { test, expect } from './fixtures.js';

test('timeline mode keeps camera rotation and zoom responsive', async ({ page, loadViewerPage }) => {
  test.setTimeout(120_000);
  await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

  await page.evaluate(() => {
    window.viewerApi.startMDContinuous({ steps: 220, temperature: 1400 });
  });

  await page.waitForFunction((target) => {
    const stats = window.viewerApi.timeline.bufferStats();
    return stats.size >= target;
  }, 50, { timeout: 45_000 });

  await page.evaluate(() => window.viewerApi.stopSimulation());

  const hitbox = page.locator('[data-testid="timeline-hitbox"]');
  await hitbox.hover();

  const slider = page.locator('[data-testid="timeline-slider"]');
  await slider.evaluate((node) => {
    node.value = '-12';
    node.dispatchEvent(new Event('input', { bubbles: true }));
    node.dispatchEvent(new Event('change', { bubbles: true }));
  });

  await page.waitForFunction(() => window.viewerApi.timeline.getState().active, null, { timeout: 12_000 });

  const dragMeta = await page.evaluate(() => {
    const api = window.viewerApi;
    const canvas = api?.scene?.getEngine?.()?.getRenderingCanvas?.();
    if (!canvas) return null;
    const rect = canvas.getBoundingClientRect();
    return {
      centerX: rect.left + rect.width / 2,
      centerY: rect.top + rect.height / 2,
      beforeAlpha: api.camera.alpha,
      beforeBeta: api.camera.beta,
    };
  });
  expect(dragMeta).toBeTruthy();

  const dragDx = 160;
  const dragDy = 90;
  await page.mouse.move(dragMeta.centerX, dragMeta.centerY);
  await page.mouse.down();
  await page.mouse.move(dragMeta.centerX + dragDx, dragMeta.centerY + dragDy, { steps: 15 });
  await page.mouse.up();
  await page.waitForTimeout(50);

  const rotateAfter = await page.evaluate(() => {
    const api = window.viewerApi;
    return { alpha: api.camera.alpha, beta: api.camera.beta };
  });
  expect(Math.abs(rotateAfter.alpha - dragMeta.beforeAlpha)).toBeGreaterThan(0.001);

  const zoomBefore = await page.evaluate(() => window.viewerApi.camera.radius);
  await page.mouse.move(dragMeta.centerX, dragMeta.centerY);
  await page.mouse.wheel(0, -400);
  await page.waitForTimeout(50);
  const zoomAfter = await page.evaluate(() => window.viewerApi.camera.radius);
  expect(Math.abs(zoomAfter - zoomBefore)).toBeGreaterThan(0.001);

  const dragAttempt = await page.evaluate(() => window.viewerApi.manipulation.beginDrag(
    () => ({ x: 0, y: 0, z: 0 }),
    { planePoint: { x: 0, y: 0, z: 0 }, planeNormal: { x: 1, y: 0, z: 0 } },
  ));
  expect(dragAttempt).toBe(false);

  await page.click('[data-testid="timeline-live"]');
  await page.waitForFunction(() => !window.viewerApi.timeline.getState().active, null, { timeout: 20_000 });
});
