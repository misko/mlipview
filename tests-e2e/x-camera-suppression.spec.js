import { test, expect } from './fixtures.js';

test('x-camera detaches during drag interactions', async ({ page, loadViewerPage }) => {
  await loadViewerPage({
    query: { mol: 'molecules/water.xyz', autoMD: 0 },
    disableAutoMd: true,
  });

  await page.evaluate(() => {
    window.viewerApi.setTestAutoSelectFallback(true);
    window.viewerApi.debugCameraControls({ reset: true });
  });

  const canvas = page.locator('#viewer');
  await expect(canvas).toBeVisible();
  const box = await canvas.boundingBox();
  expect(box).not.toBeNull();
  const { x, y, width, height } = box;
  const centerX = x + width / 2;
  const centerY = y + height / 2;

  await page.mouse.move(centerX, centerY);
  await page.mouse.down();
  await page.mouse.move(centerX + 20, centerY + 20);
  await page.mouse.up();
  await page.waitForTimeout(120);

  const stats = await page.evaluate(() => window.viewerApi.debugCameraControls());
  expect(stats.detachCalls).toBeGreaterThan(0);
  expect(stats.attachCalls).toBeGreaterThan(0);
});
