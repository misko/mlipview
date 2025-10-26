import { test, expect } from './fixtures.js';

test.describe('UI VR integration and canvas sizing', () => {
  test('does not render Babylon default XR button', async ({ page, loadViewerPage }) => {
    test.setTimeout(45000);
    await loadViewerPage({ query: { autoMD: 0 }, testMode: false });
    const vrIcon = page.locator('.babylonVRicon');
    await expect(vrIcon).toHaveCount(0, { timeout: 4000 });
  });

  test('canvas resize keeps render resolution in sync with viewport', async ({
    page,
    loadViewerPage,
  }) => {
    test.setTimeout(45000);
    await loadViewerPage({ query: { autoMD: 0 }, testMode: false });
    await page.waitForFunction(
      () => typeof window !== 'undefined' && !!window.viewerApi && !!window.viewerApi.scene,
    );

    const sample = async () =>
      await page.evaluate(() => {
        const canvas = document.getElementById('viewer');
        if (!canvas) return null;
        return {
          width: canvas.width,
          height: canvas.height,
          clientWidth: canvas.clientWidth,
          clientHeight: canvas.clientHeight,
          dpr: window.devicePixelRatio || 1,
        };
      });

    await page.setViewportSize({ width: 1200, height: 800 });
    await page.waitForTimeout(200);
    const before = await sample();
    expect(before).not.toBeNull();
    if (!before) throw new Error('viewer canvas not ready before resize');

    await page.setViewportSize({ width: 900, height: 600 });
    await page.waitForFunction(() => {
      const canvas = document.getElementById('viewer');
      if (!canvas) return false;
      const dpr = window.devicePixelRatio || 1;
      const expectedW = Math.round(canvas.clientWidth * dpr);
      const expectedH = Math.round(canvas.clientHeight * dpr);
      return Math.abs(canvas.width - expectedW) <= 2 && Math.abs(canvas.height - expectedH) <= 2;
    });
    const after = await sample();
    expect(after).not.toBeNull();
    if (!after) throw new Error('viewer canvas not ready after resize');

    const first = before;
    const second = after;
    expect(second.width).not.toBe(first.width);
    expect(second.height).not.toBe(first.height);
    const tolerance = 3;
    expect(Math.abs(second.width - second.clientWidth * second.dpr)).toBeLessThanOrEqual(tolerance);
    expect(Math.abs(second.height - second.clientHeight * second.dpr)).toBeLessThanOrEqual(tolerance);
  });
});
