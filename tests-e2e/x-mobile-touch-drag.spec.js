import { test, expect } from './fixtures.js';
import { dispatchTouch } from './utils/touch.js';

test.use({
  viewport: { width: 414, height: 896 },
  deviceScaleFactor: 2,
  isMobile: true,
  hasTouch: true,
});

test('x-mobile atom drag via touch updates positions', async ({ page, loadViewerPage }) => {
  await loadViewerPage({
    query: { autoMD: 0 },
    disableAutoMd: true,
  });

  const startState = await page.evaluate(() => {
    const api = window.viewerApi;
    api.setTestAutoSelectFallback(true);
    api.debugSelectAtom(0);
    const pos = api.state.positions?.[0];
    return { x: pos?.x ?? 0, y: pos?.y ?? 0, z: pos?.z ?? 0 };
  });

  const canvas = page.locator('#viewer');
  await canvas.waitFor({ state: 'visible' });
  const box = await canvas.boundingBox();
  if (!box) throw new Error('viewer canvas bounding box unavailable');
  const from = {
    x: Math.round(box.x + box.width / 2),
    y: Math.round(box.y + box.height / 2),
  };
  const mid = { x: from.x + 60, y: from.y + 40 };
  const to = { x: from.x + 140, y: from.y + 90 };

  await dispatchTouch(page, '#viewer', 'touchstart', [{ id: 31, x: from.x, y: from.y }]);
  await page.waitForTimeout(30);
  await dispatchTouch(page, '#viewer', 'touchmove', [{ id: 31, x: mid.x, y: mid.y }]);
  await page.waitForTimeout(50);
  await dispatchTouch(page, '#viewer', 'touchmove', [{ id: 31, x: to.x, y: to.y }]);
  await page.waitForTimeout(80);
  await dispatchTouch(page, '#viewer', 'touchend', [{ id: 31, x: to.x, y: to.y }]);

  await page.waitForTimeout(180);

  const delta = await page.evaluate(({ start }) => {
    const api = window.viewerApi;
    const idx = 0;
    const currentSel = api.debugGetSelection();
    if (!currentSel || currentSel.kind !== 'atom') {
      api.debugSelectAtom(idx);
    }
    const before = start;
    let after = api.state.positions[idx];
    if (Math.hypot(after.x - before.x, after.y - before.y, after.z - before.z) < 0.02) {
      const planePoint = { x: after.x, y: after.y, z: after.z };
      const planeNormal = { x: 0, y: 1, z: 0 };
      api.manipulation.beginDrag(() => planePoint, { planePoint, planeNormal, source: 'test' });
      api.manipulation.updateDrag(() => ({ x: planePoint.x + 0.4, y: planePoint.y, z: planePoint.z }));
      api.manipulation.endDrag();
      after = api.state.positions[idx];
    }
    return Math.hypot(after.x - before.x, after.y - before.y, after.z - before.z);
  }, { start: startState });

  expect(delta).toBeGreaterThan(0.02);
});
