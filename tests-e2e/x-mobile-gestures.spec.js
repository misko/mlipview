import { test, expect } from './fixtures.js';
import { dispatchTouch } from './utils/touch.js';

test.use({
  viewport: { width: 414, height: 896 },
  deviceScaleFactor: 2,
  isMobile: true,
  hasTouch: true,
});

async function getCanvasCenter(page) {
  const canvas = page.locator('#viewer');
  await canvas.waitFor({ state: 'visible' });
  const box = await canvas.boundingBox();
  if (!box) throw new Error('viewer canvas bounding box unavailable');
  return {
    x: Math.round(box.x + box.width / 2),
    y: Math.round(box.y + box.height / 2),
  };
}

test('x-mobile single-finger rotate does not zoom', async ({ page, loadViewerPage }) => {
  await loadViewerPage({
    query: { autoMD: 0 },
    disableAutoMd: true,
  });

  const before = await page.evaluate(() => ({
    alpha: window.viewerApi.camera?.alpha ?? 0,
    beta: window.viewerApi.camera?.beta ?? 0,
    radius: window.viewerApi.camera?.radius ?? 0,
  }));

  const center = await getCanvasCenter(page);

  await dispatchTouch(page, '#viewer', 'touchstart', [{ id: 11, x: center.x, y: center.y }]);
  await dispatchTouch(page, '#viewer', 'touchmove', [
    { id: 11, x: center.x + 120, y: center.y + 20 },
  ]);
  await dispatchTouch(page, '#viewer', 'touchend', [
    { id: 11, x: center.x + 120, y: center.y + 20 },
  ]);

  await page.waitForTimeout(80);

  const after = await page.evaluate(() => ({
    alpha: window.viewerApi.camera?.alpha ?? 0,
    beta: window.viewerApi.camera?.beta ?? 0,
    radius: window.viewerApi.camera?.radius ?? 0,
  }));

  expect(Math.abs(after.alpha - before.alpha)).toBeGreaterThan(0.1);
  expect(Math.abs(after.radius - before.radius)).toBeLessThan(0.05);
});

test('x-mobile pinch gesture changes camera zoom', async ({ page, loadViewerPage }) => {
  await loadViewerPage({
    query: { autoMD: 0 },
    disableAutoMd: true,
  });

  const beforeRadius = await page.evaluate(() => window.viewerApi.camera?.radius ?? 0);
  const center = await getCanvasCenter(page);

  const start = [
    { id: 21, x: center.x - 80, y: center.y },
    { id: 22, x: center.x + 80, y: center.y },
  ];
  const move = [
    { id: 21, x: center.x - 140, y: center.y },
    { id: 22, x: center.x + 140, y: center.y },
  ];

  await dispatchTouch(page, '#viewer', 'touchstart', start);
  await dispatchTouch(page, '#viewer', 'touchmove', move);
  await dispatchTouch(page, '#viewer', 'touchend', move);

  await page.waitForTimeout(80);

  const afterRadius = await page.evaluate(() => window.viewerApi.camera?.radius ?? 0);
  expect(Math.abs(afterRadius - beforeRadius)).toBeGreaterThan(0.3);
});
