import { test, expect } from './fixtures.js';

test.describe('timeline slider selection', () => {
  test('single slider change selects requested frame', async ({ page, loadViewerPage }) => {
    test.setTimeout(150_000);
    await loadViewerPage({ query: { mol: 'molecules/water.xyz' } });

    await page.evaluate(() => {
      window.viewerApi.startMDContinuous({ steps: 220, temperature: 1600 });
    });

    await page.waitForFunction(
      (target) => {
        const stats = window.viewerApi.timeline.bufferStats();
        return stats.size >= target;
      },
      40,
      { timeout: 45_000 },
    );

    const hitbox = page.locator('[data-testid="timeline-hitbox"]');
    await hitbox.hover();
    const slider = page.locator('[data-testid="timeline-slider"]');
    await slider.waitFor({ state: 'visible' });

    await page.waitForFunction(() => {
      const sliderEl = document.querySelector('[data-testid="timeline-slider"]');
      return sliderEl && Number(sliderEl.min) < -1;
    }, null, { timeout: 20_000 });

    await page.evaluate(() => {
      const sliderEl = document.querySelector('[data-testid="timeline-slider"]');
      if (!sliderEl) return;
      window.__pdCount = 0;
      window.__mdCount = 0;
      window.__changeCount = 0;
      sliderEl.addEventListener('pointerdown', () => { window.__pdCount++; });
      sliderEl.addEventListener('mousedown', () => { window.__mdCount++; });
      sliderEl.addEventListener('change', () => { window.__changeCount++; });
    });

    const box = await slider.boundingBox();
    expect(box).toBeTruthy();
    const clickX = box.x + box.width * 0.35;
    const clickY = box.y + box.height / 2;
    await page.mouse.click(clickX, clickY);

    const counts = await page.evaluate(() => ({
      pd: window.__pdCount || 0,
      md: window.__mdCount || 0,
      change: window.__changeCount || 0,
      state: window.viewerApi.timeline.getState(),
    }));
    console.log('counts after click', counts);

    await page.waitForFunction(
      () => {
        const state = window.viewerApi.timeline.getState();
        return state.active && state.mode === 'paused';
      },
      { timeout: 10_000 },
    );

    const state = await page.evaluate(() => window.viewerApi.timeline.getState());
    expect(state.offset).toBeLessThan(-1);
    expect(state.mode).toBe('paused');
  });
});
