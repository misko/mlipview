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

  const rotateResult = await page.evaluate(() => {
    const api = window.viewerApi;
    const canvas = api?.scene?.getEngine?.()?.getRenderingCanvas?.();
    if (!canvas) return null;
    const rect = canvas.getBoundingClientRect();
    const startX = rect.left + rect.width / 2;
    const startY = rect.top + rect.height / 2;
    const endX = startX + Math.max(32, rect.width * 0.12);
    const endY = startY + Math.max(12, rect.height * 0.06);
    const before = { alpha: api.camera.alpha, beta: api.camera.beta };
    const base = { bubbles: true, cancelable: true, pointerId: 42, pointerType: 'mouse' };
    const dispatch = (type, x, y, extras = {}) => {
      const evt = new PointerEvent(type, { ...base, clientX: x, clientY: y, ...extras });
      canvas.dispatchEvent(evt);
    };
    dispatch('pointerdown', startX, startY, { buttons: 1 });
    dispatch('pointermove', endX, endY, { buttons: 1 });
    dispatch('pointerup', endX, endY, { buttons: 0 });
    return { before, after: { alpha: api.camera.alpha, beta: api.camera.beta } };
  });

  expect(rotateResult).toBeTruthy();
  expect(Math.abs(rotateResult.after.alpha - rotateResult.before.alpha)).toBeGreaterThan(0.001);

  const zoomResult = await page.evaluate(() => {
    const api = window.viewerApi;
    const canvas = api?.scene?.getEngine?.()?.getRenderingCanvas?.();
    if (!canvas) return null;
    const rect = canvas.getBoundingClientRect();
    const centerX = rect.left + rect.width / 2;
    const centerY = rect.top + rect.height / 2;
    const before = api.camera.radius;
    const wheel = new WheelEvent('wheel', {
      bubbles: true,
      cancelable: true,
      deltaY: -240,
      clientX: centerX,
      clientY: centerY,
    });
    canvas.dispatchEvent(wheel);
    return { before, after: api.camera.radius };
  });

  expect(zoomResult).toBeTruthy();
  expect(Math.abs(zoomResult.after - zoomResult.before)).toBeGreaterThan(0.001);

  const dragAttempt = await page.evaluate(() => window.viewerApi.manipulation.beginDrag(
    () => ({ x: 0, y: 0, z: 0 }),
    { planePoint: { x: 0, y: 0, z: 0 }, planeNormal: { x: 1, y: 0, z: 0 } },
  ));
  expect(dragAttempt).toBe(false);

  await page.click('[data-testid="timeline-live"]');
  await page.waitForFunction(() => !window.viewerApi.timeline.getState().active, null, { timeout: 20_000 });
});
