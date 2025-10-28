import { test, expect } from './fixtures.js';
import { dispatchTouch } from './utils/touch.js';

test.use({
  viewport: { width: 414, height: 896 },
  deviceScaleFactor: 2,
  isMobile: true,
  hasTouch: true,
});

test('x-mobile camera detaches and reattaches around touch gestures', async ({
  page,
  loadViewerPage,
}) => {
  await loadViewerPage({
    query: { autoMD: 0 },
    disableAutoMd: true,
    configOverrides: { dragThrottleMs: 16 },
  });

  await page.evaluate(() => {
    window.viewerApi.debugCameraControls({ reset: true });
  });

  const canvas = page.locator('#viewer');
  await canvas.waitFor({ state: 'visible' });
  const box = await canvas.boundingBox();
  if (!box) throw new Error('viewer canvas bounding box unavailable');
  const center = {
    x: Math.round(box.x + box.width / 2),
    y: Math.round(box.y + box.height / 2),
  };

  await dispatchTouch(page, '#viewer', 'touchstart', [{ id: 1, x: center.x, y: center.y }]);
  await dispatchTouch(page, '#viewer', 'touchmove', [{ id: 1, x: center.x + 40, y: center.y + 32 }]);
  await dispatchTouch(page, '#viewer', 'touchend', [{ id: 1, x: center.x + 40, y: center.y + 32 }]);

  await page.waitForTimeout(50);

  const stats = await page.evaluate(() => window.viewerApi.debugCameraControls());
  expect(stats.detachCalls).toBeGreaterThan(0);
  expect(stats.attachCalls).toBeGreaterThan(0);
});
