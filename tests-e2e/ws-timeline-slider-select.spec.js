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

    const targetOffset = await page.evaluate(() => {
      const offsets = window.viewerApi.timeline.getOffsets();
      if (!offsets || !offsets.length) return null;
      const candidates = offsets.filter((o) => o < -1);
      return candidates.length ? candidates[Math.floor(candidates.length / 2)] : null;
    });
    expect(targetOffset).not.toBeNull();

    await page.evaluate(() => {
      window.__changeCount = 0;
      const sliderEl = document.querySelector('[data-testid="timeline-slider"]');
      sliderEl?.addEventListener('change', () => { window.__changeCount++; }, { once: false });
    });

    await slider.evaluate((node, value) => {
      if (!node) return;
      node.value = String(value);
      node.dispatchEvent(new Event('input', { bubbles: true }));
      node.dispatchEvent(new Event('change', { bubbles: true }));
    }, targetOffset);

    await page.waitForFunction((value) => {
      const state = window.viewerApi.timeline.getState();
      return state.offset === value && state.mode !== 'live';
    }, targetOffset, { timeout: 10_000 });

    const result = await page.evaluate(() => ({
      state: window.viewerApi.timeline.getState(),
      change: window.__changeCount || 0,
    }));
    expect(result.state.offset).toBe(targetOffset);
    expect(result.state.mode).toBe('paused');
    expect(result.change).toBeGreaterThanOrEqual(1);
  });
});
